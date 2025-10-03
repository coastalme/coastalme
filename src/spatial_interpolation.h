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

   // Get the k-d tree (for sharing between interpolators)
   KDTree const* GetKDTree() const { return m_kdtree; }

   // Get the point cloud (for sharing between interpolators)
   PointCloud const& GetPointCloud() const { return m_cloud; }

   private:
   PointCloud m_cloud;
   std::vector<double> m_values;
   KDTree* m_kdtree;
   int m_k_neighbors;
   double m_power;
   bool m_owns_kdtree;

   static constexpr double EPSILON = 1e-10;

   // Private constructor for sharing k-d tree
   friend class DualSpatialInterpolator;
   SpatialInterpolator(PointCloud const& cloud,
                       KDTree* kdtree,
                       std::vector<double> const& values,
                       int k_neighbors,
                       double power);
};

// Optimized dual interpolator for X and Y values sharing same spatial points
class DualSpatialInterpolator
{
   public:
   DualSpatialInterpolator(std::vector<Point2D> const& points,
                           std::vector<double> const& values_x,
                           std::vector<double> const& values_y,
                           int k_neighbors = 12,
                           double power = 2.0);

   ~DualSpatialInterpolator();

   // Interpolate both X and Y at multiple points (parallel optimized)
   void Interpolate(std::vector<Point2D> const& query_points,
                    std::vector<double>& results_x,
                    std::vector<double>& results_y) const;

   private:
   PointCloud m_cloud;
   std::vector<double> m_values_x;
   std::vector<double> m_values_y;
   KDTree* m_kdtree;
   int m_k_neighbors;
   double m_power;

   static constexpr double EPSILON = 1e-10;

   // Helper for single point interpolation
   void InterpolatePoint(double x, double y, double& result_x, double& result_y,
                         std::vector<unsigned int>& indices,
                         std::vector<double>& sq_dists) const;
};

#endif      // SPATIAL_INTERPOLATION_H
