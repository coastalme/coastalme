#ifndef SPATIAL_INTERPOLATION_H
#define SPATIAL_INTERPOLATION_H

#include <vector>
#include <cmath>
#include <algorithm>
#include "nanoflann.hpp"

// Simple 2D point structure
struct Point2D
{
   double x, y;
   Point2D(double x_ = 0, double y_ = 0) : x(x_), y(y_)
   {
   }
};

// Point cloud adaptor for nanoflann
struct PointCloud
{
   std::vector<Point2D> pts;

   inline size_t kdtree_get_point_count() const
   {
      return pts.size();
   }

   inline double kdtree_get_pt(size_t const idx, size_t const dim) const
   {
      return (dim == 0) ? pts[idx].x : pts[idx].y;
   }

   template<class BBOX>
   bool kdtree_get_bbox(BBOX&) const
   {
      return false;
   }
};

// Type definitions for nanoflann
using KDTree = nanoflann::KDTreeSingleIndexAdaptor<
   nanoflann::L2_Simple_Adaptor<double, PointCloud>,
   PointCloud,
   2      // 2D
   >;

// Main interpolator class
class SpatialInterpolator
{
   public:
   SpatialInterpolator(std::vector<Point2D> const& points,
                       std::vector<double> const& values,
                       int k_neighbors = 12,
                       double power = 2.0);

   ~SpatialInterpolator();

   // Interpolate at a single point
   double Interpolate(double x, double y) const;

   // Interpolate at multiple points (more efficient)
   void Interpolate(std::vector<Point2D> const& query_points,
                    std::vector<double>& results) const;

   private:
   PointCloud m_cloud;
   std::vector<double> m_values;
   KDTree* m_kdtree;
   int m_k_neighbors;
   double m_power;

   static constexpr double EPSILON = 1e-10;
};

#endif      // SPATIAL_INTERPOLATION_H
