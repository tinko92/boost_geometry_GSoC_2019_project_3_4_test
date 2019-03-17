#include <random>
#include <vector>
#include <algorithm>
#include <iterator>

#include <boost/geometry/geometry.hpp>

using namespace boost::geometry;

template
<
    typename OutPoint,
    typename InPoint,
    typename Generator
>
OutPoint sample_triangle_uniform(const InPoint& p1, const InPoint& p2, const InPoint& p3, Generator& gen)
{
    std::uniform_real_distribution<> dis(0, 1);
    double r1 = dis(gen);
    double r2 = dis(gen);
    return make<OutPoint>(
        (1-std::sqrt(r1))*get<0>(p1)+
        (std::sqrt(r1)*(1-r2))*get<0>(p2)+
        (r2*std::sqrt(r1))*get<0>(p3),
        (1-std::sqrt(r1))*get<1>(p1)+
        (std::sqrt(r1)*(1-r2))*get<1>(p2)+
        (r2*std::sqrt(r1))*get<1>(p3));
}

template
<
    typename Polygon,
    typename Collection
>
void sample_polygon_uniform(const Polygon& p, std::size_t num, Collection& out)
{
    typedef typename point_type<Polygon>::type point_type;
    typedef typename area_result<Polygon>::type area_type;
    typedef typename boost::range_value<Collection>::type out_point_type;
    typename ring_type<const Polygon>::type& ring = exterior_ring(p);
    const std::size_t length = boost::size(ring);
    std::vector<typename area_result<Polygon>::type> cumulative_area;
    cumulative_area.reserve(length-2);
    cumulative_area.push_back(0);
    area_type cum_sum = 0;
    for(std::size_t i = 1;i<length-2;++i) {
        cum_sum += triangle_double_area<area_type>(*boost::begin(ring), *(boost::begin(ring)+i), *(boost::begin(ring)+i+1));
        cumulative_area.push_back(cum_sum);
    }

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0, 1);
    for(std::size_t i = 0; i<num; ++i) {
        area_type r = boost::numeric_cast<area_type>(dis(gen))*cum_sum;
        std::size_t t = std::distance(cumulative_area.begin(),
           std::lower_bound(cumulative_area.begin(),cumulative_area.end(), r));
        *std::back_inserter(out) = sample_triangle_uniform<out_point_type>(*boost::begin(ring), *(boost::begin(ring)+t), *(boost::begin(ring)+t+1),gen);
    }
}
