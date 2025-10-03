#include "spatial_interpolation.h"
#include <stdexcept>

#ifdef _OPENMP
#include <omp.h>
#endif

SpatialInterpolator::SpatialInterpolator(std::vector<Point2D> const& points,
                                         std::vector<double> const& values,
                                         int k_neighbors,
                                         double power)
   : m_values(values), m_kdtree(nullptr), m_k_neighbors(k_neighbors), m_power(power), m_owns_kdtree(true)
{
   if (points.size() != values.size())
      throw std::invalid_argument("Points and values must have same size");

   if (points.empty())
      throw std::invalid_argument("Cannot create interpolator with empty data");

   // Copy points
   m_cloud.pts = points;

   // Build k-d tree
   m_kdtree = new KDTree(2, m_cloud, {10 /* max leaf size */});
   m_kdtree->buildIndex();
}

SpatialInterpolator::SpatialInterpolator(PointCloud const& cloud,
                                         KDTree* kdtree,
                                         std::vector<double> const& values,
                                         int k_neighbors,
                                         double power)
   : m_cloud(cloud), m_values(values), m_kdtree(kdtree), m_k_neighbors(k_neighbors), m_power(power), m_owns_kdtree(false)
{
}

SpatialInterpolator::~SpatialInterpolator()
{
   if (m_owns_kdtree)
      delete m_kdtree;
}

double SpatialInterpolator::Interpolate(double x, double y) const
{
   // Prepare query
   double const query_pt[2] = {x, y};

   // Find k nearest neighbors
   size_t const k = std::min((size_t) m_k_neighbors, m_cloud.pts.size());
   std::vector<unsigned int> indices(k);
   std::vector<double> sq_dists(k);

   unsigned int num_found = m_kdtree->knnSearch(query_pt, k,
                                                 indices.data(),
                                                 sq_dists.data());

   if (num_found == 0)
      throw std::runtime_error("knnSearch found no neighbors");

   // Check if we're exactly on a data point
   if (sq_dists[0] < EPSILON)
      return m_values[indices[0]];

   // Inverse Distance Weighting (IDW) - optimized for power=2.0
   double sum_weights = 0.0;
   double sum_weighted_values = 0.0;

   if (m_power == 2.0)
   {
      // Optimized path: use squared distances directly, avoid sqrt and pow
      for (size_t i = 0; i < num_found; i++)
      {
         double weight = 1.0 / sq_dists[i];  // 1/dist^2 = 1/sq_dist
         sum_weights += weight;
         sum_weighted_values += weight * m_values[indices[i]];
      }
   }
   else
   {
      // General case
      for (size_t i = 0; i < num_found; i++)
      {
         double dist = std::sqrt(sq_dists[i]);
         double weight = 1.0 / std::pow(dist, m_power);
         sum_weights += weight;
         sum_weighted_values += weight * m_values[indices[i]];
      }
   }

   return sum_weighted_values / sum_weights;
}

void SpatialInterpolator::Interpolate(std::vector<Point2D> const& query_points,
                                      std::vector<double>& results) const
{
   results.resize(query_points.size());
   size_t const k = std::min((size_t) m_k_neighbors, m_cloud.pts.size());

#ifdef _OPENMP
   #pragma omp parallel
   {
      // Thread-local storage to avoid allocation overhead
      std::vector<unsigned int> indices(k);
      std::vector<double> sq_dists(k);

      #pragma omp for schedule(guided, 128) nowait
      for (size_t i = 0; i < query_points.size(); i++)
      {
         double const query_pt[2] = {query_points[i].x, query_points[i].y};

         unsigned int num_found = m_kdtree->knnSearch(query_pt, k,
                                                       indices.data(),
                                                       sq_dists.data());

         if (num_found == 0)
         {
            results[i] = 0.0;  // or m_dMissingValue
            continue;
         }

         // Check if we're exactly on a data point
         if (sq_dists[0] < EPSILON)
         {
            results[i] = m_values[indices[0]];
            continue;
         }

         // IDW calculation - optimized for power=2.0
         double sum_weights = 0.0;
         double sum_weighted_values = 0.0;

         if (m_power == 2.0)
         {
            for (size_t j = 0; j < num_found; j++)
            {
               double weight = 1.0 / sq_dists[j];
               sum_weights += weight;
               sum_weighted_values += weight * m_values[indices[j]];
            }
         }
         else
         {
            for (size_t j = 0; j < num_found; j++)
            {
               double dist = std::sqrt(sq_dists[j]);
               double weight = 1.0 / std::pow(dist, m_power);
               sum_weights += weight;
               sum_weighted_values += weight * m_values[indices[j]];
            }
         }

         results[i] = sum_weighted_values / sum_weights;
      }
   }
#else
   // Serial fallback - still optimized with reused buffers
   std::vector<unsigned int> indices(k);
   std::vector<double> sq_dists(k);

   for (size_t i = 0; i < query_points.size(); i++)
   {
      results[i] = Interpolate(query_points[i].x, query_points[i].y);
   }
#endif
}

// DualSpatialInterpolator implementation
DualSpatialInterpolator::DualSpatialInterpolator(std::vector<Point2D> const& points,
                                                 std::vector<double> const& values_x,
                                                 std::vector<double> const& values_y,
                                                 int k_neighbors,
                                                 double power)
   : m_values_x(values_x), m_values_y(values_y), m_kdtree(nullptr), m_k_neighbors(k_neighbors), m_power(power)
{
   if (points.size() != values_x.size() || points.size() != values_y.size())
      throw std::invalid_argument("Points and values must have same size");

   if (points.empty())
      throw std::invalid_argument("Cannot create interpolator with empty data");

   // Copy points
   m_cloud.pts = points;

   // Build k-d tree (shared for both X and Y)
   m_kdtree = new KDTree(2, m_cloud, {10 /* max leaf size */});
   m_kdtree->buildIndex();
}

DualSpatialInterpolator::~DualSpatialInterpolator()
{
   delete m_kdtree;
}

void DualSpatialInterpolator::InterpolatePoint(double x, double y,
                                               double& result_x, double& result_y,
                                               std::vector<unsigned int>& indices,
                                               std::vector<double>& sq_dists) const
{
   double const query_pt[2] = {x, y};
   size_t const k = std::min((size_t) m_k_neighbors, m_cloud.pts.size());

   unsigned int num_found = m_kdtree->knnSearch(query_pt, k,
                                                 indices.data(),
                                                 sq_dists.data());

   if (num_found == 0)
   {
      result_x = result_y = 0.0;
      return;
   }

   // Check if we're exactly on a data point
   if (sq_dists[0] < EPSILON)
   {
      result_x = m_values_x[indices[0]];
      result_y = m_values_y[indices[0]];
      return;
   }

   // IDW for both X and Y simultaneously - optimized for power=2.0
   double sum_weights = 0.0;
   double sum_weighted_x = 0.0;
   double sum_weighted_y = 0.0;

   if (m_power == 2.0)
   {
      for (size_t i = 0; i < num_found; i++)
      {
         double weight = 1.0 / sq_dists[i];
         sum_weights += weight;
         sum_weighted_x += weight * m_values_x[indices[i]];
         sum_weighted_y += weight * m_values_y[indices[i]];
      }
   }
   else
   {
      for (size_t i = 0; i < num_found; i++)
      {
         double dist = std::sqrt(sq_dists[i]);
         double weight = 1.0 / std::pow(dist, m_power);
         sum_weights += weight;
         sum_weighted_x += weight * m_values_x[indices[i]];
         sum_weighted_y += weight * m_values_y[indices[i]];
      }
   }

   result_x = sum_weighted_x / sum_weights;
   result_y = sum_weighted_y / sum_weights;
}

void DualSpatialInterpolator::Interpolate(std::vector<Point2D> const& query_points,
                                          std::vector<double>& results_x,
                                          std::vector<double>& results_y) const
{
   size_t const n = query_points.size();
   results_x.resize(n);
   results_y.resize(n);

   // Early exit for empty query
   if (n == 0) return;

   size_t const k = std::min((size_t) m_k_neighbors, m_cloud.pts.size());

#ifdef _OPENMP
   // Use guided scheduling to reduce overhead
   #pragma omp parallel
   {
      // Thread-local storage
      std::vector<unsigned int> indices(k);
      std::vector<double> sq_dists(k);

      #pragma omp for schedule(guided, 128) nowait
      for (size_t i = 0; i < query_points.size(); i++)
      {
         InterpolatePoint(query_points[i].x, query_points[i].y,
                         results_x[i], results_y[i],
                         indices, sq_dists);
      }
   }
#else
   // Serial fallback
   std::vector<unsigned int> indices(k);
   std::vector<double> sq_dists(k);

   for (size_t i = 0; i < query_points.size(); i++)
   {
      InterpolatePoint(query_points[i].x, query_points[i].y,
                      results_x[i], results_y[i],
                      indices, sq_dists);
   }
#endif
}
