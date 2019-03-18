#include <iostream>
#include <fstream>
#include <random>
#include <chrono>

#include "shull.hpp"
#include "convex-polygon-sampler.hpp"

namespace bg = boost::geometry;
typedef bg::model::point<double, 2, bg::cs::cartesian> point_t;
typedef bg::model::segment<point_t> segment_t;
typedef bg::model::box<point_t> box_t;
typedef bg::model::polygon<point_t, false> polygon_t;
typedef bg::model::multi_point<point_t> mpoint_t;
typedef triangle<point_t> Triangle;

int main() {
	mpoint_t mpt2;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> dis(200, 100);

    //larger instance for performance measurement
    std::size_t num_points = 100'000;

    for(std::size_t i = 0; i<num_points; ++i)
        bg::append(mpt2, point_t(dis(gen), dis(gen)));

    auto start = std::chrono::system_clock::now();

	shull(mpt2);

    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Time for Delaunay-Triangulation of " << num_points << " random points: " << elapsed.count() << "ms\n";

    //small instance for drawing
    num_points = 1'000;
    mpt2 = {};

    for(std::size_t i = 0; i<num_points; ++i)
        bg::append(mpt2, point_t(dis(gen), dis(gen)));

    auto [triangulation,hull] = shull(mpt2);

    std::ofstream svg("triangulation.svg");
    bg::svg_mapper<point_t> mapper(svg, 800, 800);
    for(const auto& t : triangulation) {
        mapper.add(segment_t(t.p1,t.p2));
        mapper.add(segment_t(t.p2,t.p3));
        mapper.add(segment_t(t.p3,t.p1));
    }
    for(const auto& t : triangulation) {
        mapper.map(segment_t(t.p1,t.p2), "stroke:rgb(153,204,0);stroke-width:0.1",1);
        mapper.map(segment_t(t.p2,t.p3), "stroke:rgb(153,204,0);stroke-width:0.1",1);
        mapper.map(segment_t(t.p3,t.p1), "stroke:rgb(153,204,0);stroke-width:0.1",1);
    }
    for(const auto& s : hull) {
        mapper.add(s);
    }
    for(const auto& s : hull) {
        mapper.map(s, "stroke:rgb(255,0,0);stroke-width:1",1);
    }

    polygon_t poly;
    bg::append(poly.outer(),hull[0].first);
    for(const auto& s : hull) {
        bg::append(poly.outer(), s.second);
    }
    std::vector<point_t> out;
    const std::size_t sample_count = 5000;

    sample_polygon_uniform(poly,sample_count,out);

    std::ofstream svg2("random_samples.svg");
    bg::svg_mapper<point_t> mapper2(svg2, 800, 800);
    mapper.add(poly);
    for(const auto& p : out) {
        mapper2.add(p);
    }
    mapper2.map(poly, "opacity:0.4;fill:none;stroke:rgb(212,0,0);stroke-width:5");
    for(const auto& p : out) {
        mapper2.map(p, "fill-opacity:0.5;fill:rgb(153,204,0);stroke:rgb(153,204,0);stroke-width:1",1);
    }

    bool inside = true;
    for(const auto& p : out) {
        inside = inside && bg::within(p,poly);
    }
    std::cout << "All points inside polygon? " << (inside?"yes":"no") << "\n";
    box_t b;
    bg::envelope(poly, b);
    box_t lower(b.min_corner(),point_t(bg::get<0>(b.max_corner()),bg::get<1>(b.max_corner())/2));
    bg::model::multi_polygon<polygon_t> lower_poly;
    bg::intersection(lower,poly,lower_poly);
    auto lower_area_ratio = bg::area(lower_poly[0].outer())/bg::area(poly);
    std::size_t count_lower = 0;
    for(const auto& p : out) {
        if(bg::within(p,lower))
            ++count_lower;
    }
    std::cout << count_lower << "/" << sample_count << " (" << boost::numeric_cast<double>(count_lower)/sample_count << ") points in lower half. Area ratio: " << lower_area_ratio << "\n";
	return 0;
}
