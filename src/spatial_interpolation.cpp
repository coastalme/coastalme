#include "spatial_interpolation.h"
#include <stdexcept>

SpatialInterpolator::SpatialInterpolator(std::vector<Point2D> const& points,
                                         std::vector<double> const& values,
                                         int k_neighbors,
                                         double power)
   : m_values(values), m_kdtree(nullptr), m_k_neighbors(k_neighbors), m_power(power)
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

SpatialInterpolator::~SpatialInterpolator()
{
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

   // Inverse Distance Weighting (IDW)
   double sum_weights = 0.0;
   double sum_weighted_values = 0.0;

   for (size_t i = 0; i < num_found; i++)
   {
      double dist = std::sqrt(sq_dists[i]);
      double weight = 1.0 / std::pow(dist, m_power);

      sum_weights += weight;
      sum_weighted_values += weight * m_values[indices[i]];
   }

   return sum_weighted_values / sum_weights;
}

void SpatialInterpolator::Interpolate(std::vector<Point2D> const& query_points,
                                      std::vector<double>& results) const
{
   results.resize(query_points.size());

   for (size_t i = 0; i < query_points.size(); i++)
   {
      results[i] = Interpolate(query_points[i].x, query_points[i].y);
   }
}
