/**
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package org.apache.lucene.spatial.geohash;

import org.apache.lucene.index.DocsEnum;
import org.apache.lucene.index.IndexReader;
import org.apache.lucene.index.Terms;
import org.apache.lucene.index.TermsEnum;
import org.apache.lucene.search.DocIdSet;
import org.apache.lucene.search.Filter;
import org.apache.lucene.spatial.DistanceUtils;
import org.apache.lucene.spatial.geometry.shape.Geometry2D;
import org.apache.lucene.spatial.geometry.shape.IntersectCase;
import org.apache.lucene.spatial.geometry.shape.Point2D;
import org.apache.lucene.spatial.geometry.shape.PointDistanceGeom;
import org.apache.lucene.util.Bits;
import org.apache.lucene.util.BytesRef;
import org.apache.lucene.util.OpenBitSet;

import java.io.IOException;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.Map;

/**
 * Performs a spatial filter against a field indexed using Geohashes. Using the hierarchical grid nature of geohashes,
 * this filter recursively traverses each precision length and uses methods on {@link Geometry2D} to efficiently know
 * that all points at a geohash fit in the shape or not to either short-circuit unnecessary traversals or to efficiently
 * load all enclosed points.
 */
public class GeoHashPrefixFilter extends Filter {

  private static final int GRIDLEN_SCAN_THRESHOLD = 4;//>= 1
  private final String fieldName;
  private final Geometry2D geoShape;
  private final GridNode.GridReferenceSystem gridReferenceSystem;

  private Map<Integer,Double> docDistanceCache;

  public GeoHashPrefixFilter(String fieldName, Geometry2D geoShape, GridNode.GridReferenceSystem gridReferenceSystem) {
    this.fieldName = fieldName;
    this.geoShape = geoShape;
    this.gridReferenceSystem = gridReferenceSystem;
  }

  @Override
  public DocIdSet getDocIdSet(IndexReader.AtomicReaderContext ctx) throws IOException {
    IndexReader reader = ctx.reader;
    OpenBitSet bits = new OpenBitSet(reader.maxDoc());
    Terms terms = reader.fields().terms(fieldName);
    if (terms == null)
      return null;
    TermsEnum termsEnum = terms.iterator();
    DocsEnum docsEnum = null;//cached for termsEnum.docs() calls
    Bits delDocs = reader.getDeletedDocs();
    BytesRef term = termsEnum.next();//the most recent term examined via termsEnum.term()

    final boolean needDistance = false;
    Point2D nd_center = null;
    double nd_radius = DistanceUtils.EARTH_MEAN_RADIUS_KM;
    double nd_maxDistance = 0;//when > 0 then this is an optimization to avoid double-haversine calculation
    if (needDistance) {
      nd_center = geoShape.centroid();
      docDistanceCache = new HashMap<Integer, Double>();
      if (geoShape instanceof PointDistanceGeom) {
        PointDistanceGeom pointDistanceGeom = (PointDistanceGeom) geoShape;
        nd_maxDistance = pointDistanceGeom.getDistance();
        nd_radius = pointDistanceGeom.getRadius();
      }
    }

    //TODO Add a precision short-circuit so that we are not accurate on the edge but we're faster.

    //TODO An array based nodes impl would be more efficient; or a stack of iterators.  LinkedList conveniently has bulk add to beginning.
    LinkedList<GridNode> nodes = new LinkedList<GridNode>(gridReferenceSystem.getSubNodes(geoShape.boundingRectangle()));
    while(!nodes.isEmpty() && term != null) {
      final GridNode node = nodes.removeFirst();
      assert node.length() > 0;
      if (!node.contains(term) && node.before(term))
        continue;//short circuit, moving >= the next indexed term
      IntersectCase intersection = geoShape.intersect(node.getRectangle());
      if (intersection == IntersectCase.OUTSIDE)
        continue;
      TermsEnum.SeekStatus seekStat = termsEnum.seek(node.getBytesRef());
      term = termsEnum.term();
      if (seekStat != TermsEnum.SeekStatus.FOUND)
        continue;
      if (intersection == IntersectCase.CONTAINS && !needDistance) {
        docsEnum = termsEnum.docs(delDocs, docsEnum);
        addDocs(docsEnum,bits, -1);
        term = termsEnum.next();//move to next term
      } else {//any other intersection
        //TODO is it worth it to optimize the shape (e.g. potentially simpler polygon)?
        //GeoShape geoShape = this.geoShape.optimize(intersection);

        //We either scan through the leaf node(s), or if there are many points then we divide & conquer.
        boolean manyPoints = node.length() < gridReferenceSystem.maxLen - GRIDLEN_SCAN_THRESHOLD;

        //TODO Try variable depth strategy:
        //IF configured to do so, we could use term.freq() as an estimate on the number of places at this depth.  OR, perhaps
        //  make estimates based on the total known term count at this level?  Or don't worry about it--use fixed depth.
//        if (manyPoints) {
//          //Make some estimations on how many points there are at this level and how few there would need to be to set
//          // manyPoints to false.
//
//          long termsThreshold = (long) estimateNumberIndexedTerms(node.length(),geoShape.getDocFreqExpenseThreshold(node));
//
//          long thisOrd = termsEnum.ord();
//          manyPoints = (termsEnum.seek(thisOrd+termsThreshold+1) != TermsEnum.SeekStatus.END
//                  && node.contains(termsEnum.term()));
//          termsEnum.seek(thisOrd);//return to last position
//        }

        if (!manyPoints) {
          //traverse all leaf terms within this node to see if they are within the geoShape, one by one.
          for(; term != null && node.contains(term); term = termsEnum.next()) {
            if (term.length < gridReferenceSystem.maxLen)//not a leaf
              continue;
            final Point2D point = gridReferenceSystem.decodeXY(term);
            //Filter those out of the shape.
            // Instead of simply calling geoShape.contains(point), which in the case of a circle calculates the distance,
            // determine if we already need the distance and if so calculate it ourselves. This is an optimization to
            // avoid calculating haversince twice.
            double distance = -1;
            if (nd_maxDistance == 0) {
              if(!geoShape.contains(point))
                continue;
            } else {
              distance = calcDistance(nd_center, point, nd_radius);
              if (distance > nd_maxDistance)
                continue;
            }
            //Calc distance if we haven't done so already
            if (needDistance) {
              if (distance == -1)
                distance = calcDistance(nd_center, point, nd_radius);
            }
            //record
            docsEnum = termsEnum.docs(delDocs, docsEnum);
            addDocs(docsEnum,bits,distance);
          }
        } else {
          //divide & conquer
          nodes.addAll(0,node.getSubNodes());//add to beginning
        }
      }
    }//node loop

    return bits;
  }

  private double calcDistance(Point2D nd_center, Point2D point, double nd_radius) {
    return DistanceUtils.haversine(Math.toRadians(nd_center.getY()), Math.toRadians(nd_center.getX()),
        Math.toRadians(point.getY()), Math.toRadians(point.getX()), nd_radius);
  }

//  double estimateNumberIndexedTerms(int nodeLen,double points) {
//    return 1000;
//    double levelProb = probabilityNumNodes[points];// [1,32)
//    if (nodeLen < geohashLength)
//      return levelProb + levelProb * estimateNumberIndexedTerms(nodeLen+1,points/levelProb);
//    return levelProb;
//  }

  private void addDocs(DocsEnum docsEnum, OpenBitSet bits, double distance) throws IOException {
    DocsEnum.BulkReadResult bulk = docsEnum.getBulkResult();
    for (; ;) {
      int nDocs = docsEnum.read();
      if (nDocs == 0) break;
      int[] docArr = bulk.docs.ints;
      int end = bulk.docs.offset + nDocs;
      for (int i = bulk.docs.offset; i < end; i++) {
        final int doc = docArr[i];
        bits.fastSet(doc);
        if (distance >= 0) {
          Double minDistance = docDistanceCache.get(doc);
          if (minDistance == null || distance < minDistance)
            docDistanceCache.put(doc,distance);
        }
      }
    }
  }

  @Override
  public String toString() {
    return "GeoFilter{fieldName='" + fieldName + '\'' + ", shape=" + geoShape + '}';
  }

  @Override
  public boolean equals(Object o) {
    if (this == o) return true;
    if (o == null || getClass() != o.getClass()) return false;

    GeoHashPrefixFilter that = (GeoHashPrefixFilter) o;

    if (fieldName != null ? !fieldName.equals(that.fieldName) : that.fieldName != null) return false;
    if (geoShape != null ? !geoShape.equals(that.geoShape) : that.geoShape != null) return false;

    return true;
  }

  @Override
  public int hashCode() {
    int result = fieldName != null ? fieldName.hashCode() : 0;
    result = 31 * result + (geoShape != null ? geoShape.hashCode() : 0);
    return result;
  }

}
