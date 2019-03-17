#include <utility>
#include <vector>
#include <stack>
#include <algorithm>

#include <boost/geometry/geometry.hpp>

using namespace boost::geometry;

template
<
    typename Point
>
struct triangle {
    Point p1,p2,p3;
    std::vector<std::size_t> adjacent_indices;
};

template
<
    typename Area,
    typename Point
>
Area triangle_double_area(const Point& p1, const Point& p2, const Point& p3) {
    return std::abs(boost::numeric_cast<Area>(get<0>(p1))*boost::numeric_cast<Area>(get<1>(p2))-boost::numeric_cast<Area>(get<0>(p1))*boost::numeric_cast<Area>(get<1>(p3))+boost::numeric_cast<Area>(get<0>(p2))*boost::numeric_cast<Area>(get<1>(p3))-boost::numeric_cast<Area>(get<0>(p2))*boost::numeric_cast<Area>(get<1>(p1))+boost::numeric_cast<Area>(get<0>(p3))*boost::numeric_cast<Area>(get<1>(p1))-boost::numeric_cast<Area>(get<0>(p3))*boost::numeric_cast<Area>(get<1>(p2)));
}

template
<
    typename Point
>
double circumcircle_diameter(const Point& p1, const Point& p2, const Point& p3)
{
    return distance(p1,p2)*distance(p1,p3)*distance(p2,p3)/triangle_double_area<double,Point>(p1,p2,p3);
}


template
<
    typename Point
>
Point circumcircle_center(const Point& p1, const Point& p2, const Point& p3)
{
    auto ax = get<0>(p1);
    auto ay = get<1>(p1);
    auto bx = get<0>(p2);
    auto by = get<1>(p2);
    auto cx = get<0>(p3);
    auto cy = get<1>(p3);
    auto d = 2*(ax*(by-cy)+bx*(cy-ay)+cx*(ay-by));
    auto x = ((ax*ax+ay*ay)*(by-cy)+(bx*bx+by*by)*(cy-ay)+(cx*cx+cy*cy)*(ay-by))/d;
    auto y = ((ax*ax+ay*ay)*(cx-bx)+(bx*bx+by*by)*(ax-cx)+(cx*cx+cy*cy)*(bx-ax))/d;
    return make<Point>(x,y);
}

template
<
    typename Segment
>
double angle(const Segment& s)
{
    return std::atan2(get<1,1>(s)-get<0,1>(s),get<1,0>(s)-get<0,0>(s));
}

template
<
    typename Point,
    typename Segment
>
bool is_left_of(const Point& p, const Segment& s)
{
    return (get<0>(p)-get<0,0>(s))*(get<1,1>(s)-get<0,1>(s))-(get<1>(p)-get<0,1>(s))*(get<1,0>(s)-get<0,0>(s))<0;
}

template
<
    typename Point
>
bool adjacent_at_edge(const triangle<Point>& t1, const triangle<Point>& t2, const Point& p1, const Point& p2)
{
    //Very naive check, not necessary with proper triangulation data structure
    return (equals(p1,t1.p1) || equals(p1,t1.p2) || equals(p1,t1.p3)) 
        && (equals(p2,t1.p1) || equals(p2,t1.p2) || equals(p2,t1.p3))
        && (equals(p1,t2.p1) || equals(p1,t2.p2) || equals(p1,t2.p3))
        && (equals(p2,t2.p1) || equals(p2,t2.p2) || equals(p2,t2.p3));
}

template
<
    typename Point
>
std::vector<std::pair<std::size_t,std::size_t>> flip_if_not_locally_delaunay(std::vector<triangle<Point>>& triangles, std::size_t ti1, std::size_t ti2)
{
    triangle<Point>& t1 = triangles[ti1];
    if(std::find(t1.adjacent_indices.begin(),t1.adjacent_indices.end(),ti2)==t1.adjacent_indices.end())
        return {};
    triangle<Point>& t2 = triangles[ti2];
    Point *p1, *p2, *p3, *p4;
    if(!equals(t1.p1,t2.p1) && !equals(t1.p1,t2.p2) && !equals(t1.p1, t2.p3)) {
        p4 = &t1.p1;
        p1 = &t1.p2;
        p3 = &t1.p3;
    }
    else if(!equals(t1.p2,t2.p1) && !equals(t1.p2,t2.p2) && !equals(t1.p2, t2.p3)) {
        p4 = &t1.p2;
        p1 = &t1.p3;
        p3 = &t1.p1;
    }
    else {
        p4 = &t1.p3;
        p1 = &t1.p1;
        p3 = &t1.p2;
    }
    if(!equals(t2.p1,*p1) && !equals(t2.p1,*p3))
        p2 = &t2.p1;
    else if(!equals(t2.p2,*p1) && !equals(t2.p2,*p3))
        p2 = &t2.p2;
    else
        p2 = &t2.p3;
    bool loc_delaunay = distance(circumcircle_center(*p1,*p2,*p3),*p4)>circumcircle_diameter(*p1,*p2,*p3)/2;
    if(loc_delaunay)
        return {};
    std::size_t t1a1=-1, t2a1=-1;
    for(const auto& i : t1.adjacent_indices) {
        if(i!=ti2)
            if(adjacent_at_edge(t1, triangles[i], *p3, *p4))
                t1a1 = i;
    }
    for(const auto& i : t2.adjacent_indices) {
        if(i!=ti1)
            if(adjacent_at_edge(t2, triangles[i], *p1, *p2))
                t2a1 = i;
    }
    Point temp;
    assign(temp,*p3);
    assign(*p3,*p2);
    if(equals(t2.p1,temp))
        assign(t2.p2,*p4);
    else if(equals(t2.p2,temp))
        assign(t2.p3,*p4);
    else
        assign(t2.p1,*p4);
    if(t1a1!=-1 && t2a1!=-1) {
        std::replace(triangles[t1a1].adjacent_indices.begin(),triangles[t1a1].adjacent_indices.end(),ti1,ti2);
        std::replace(t1.adjacent_indices.begin(),t1.adjacent_indices.end(),t1a1,t2a1);
        std::replace(triangles[t2a1].adjacent_indices.begin(),triangles[t2a1].adjacent_indices.end(),ti2,ti1);
        std::replace(t2.adjacent_indices.begin(),t2.adjacent_indices.end(),t2a1,t1a1);
    }
    else if(t1a1!=-1) {
        std::replace(triangles[t1a1].adjacent_indices.begin(),triangles[t1a1].adjacent_indices.end(),ti1,ti2);
        t1.adjacent_indices.erase(std::find(t1.adjacent_indices.begin(),t1.adjacent_indices.end(),t1a1));
        t2.adjacent_indices.push_back(t1a1);
    }
    else if(t2a1!=-1) {
        std::replace(triangles[t2a1].adjacent_indices.begin(),triangles[t2a1].adjacent_indices.end(),ti2,ti1);
        t2.adjacent_indices.erase(std::find(t2.adjacent_indices.begin(),t2.adjacent_indices.end(),t2a1));
        t1.adjacent_indices.push_back(t2a1);
    }
    std::vector<std::pair<std::size_t,std::size_t>> outer_edges;
    for(const auto& t : {ti1,ti2}) {
        for(const auto& ta : triangles[t].adjacent_indices)
            if(ta != ti1 && ta != ti2)
                outer_edges.push_back(std::pair<std::size_t,std::size_t>(t,ta));
    }
    return outer_edges;
}

template
<
    typename MultiPoint
>
std::pair<std::vector<triangle<typename boost::range_value<MultiPoint>::type>>,std::vector<model::segment<typename boost::range_value<MultiPoint>::type>>> shull(const MultiPoint& input, bool skip_delaunay = false, double eps=0.0001)
{
    typedef typename boost::range_value<MultiPoint>::type Point;
    typedef triangle<Point> Triangle;
    std::vector<Point> points;
    points.reserve(boost::size(input));
    for(const auto& in : input)
        points.push_back(in);
    //Step 2
    std::sort(points.begin()+1, points.end(), [&points](const Point& p1, const Point& p2) {
       auto d1 = comparable_distance(points[0],p1);
       auto d2 = comparable_distance(points[0],p2);
       if(d1<d2)
           return true;
       else if(d2<d1)
           return false;
       if(get<0>(p1)<get<0>(p2) || (get<0>(p1)==get<0>(p2) && get<1>(p1)<get<1>(p2) ))
           return true;
       return false;});
    for(auto it=points.end()-1;it!=points.begin();--it) //remove almost-duplicates
        if(distance(*it,*(it-1))<=eps)
            it = points.erase(it);
    //Step 4
    double min_circumcircle_diameter = circumcircle_diameter(points[0],points[1],points[2]);
    std::size_t mcd_index = 2;
    for(std::size_t i=3;i<points.size();++i) {
        auto d = circumcircle_diameter(points[0],points[1],points[i]);
        if(d<min_circumcircle_diameter) {
            mcd_index=i;
            min_circumcircle_diameter=d;
        }
        if(distance(points[0],points[i])>min_circumcircle_diameter)
            break;
    }
    std::swap(points[2],points[mcd_index]);
    Point C = circumcircle_center(points[0], points[1], points[2]);

    //Step 5
    typedef typename model::segment<Point> Segment;
    typedef typename std::pair<Segment,std::size_t> indexed_segment;
    if(!is_left_of(points[2], Segment{points[0],points[1]}))
        std::swap(points[1],points[2]);
    std::sort(points.begin()+2, points.end(), [&C](const Point& p1, const Point& p2) {
       return comparable_distance(C,p1) < comparable_distance(C,p2);
    });

    std::vector<indexed_segment> outer_hull = {indexed_segment(Segment{points[0],points[1]},0),indexed_segment(Segment{points[1],points[2]},0),indexed_segment(Segment{points[2],points[0]},0)};
    std::vector<Triangle> triangles = {Triangle{points[0],points[1],points[2],{}}};

    auto outer_hull_pred = [](const indexed_segment& s1, const indexed_segment& s2) { return angle(s1.first)<angle(s2.first); };
    //Step 6
    std::sort(outer_hull.begin(), outer_hull.end(),outer_hull_pred);

    //Step 7
    for(std::size_t i = 3; i<points.size(); ++i) {
        const Point& p = points[i];
        std::size_t end_tri = triangles.size();
        //The following could (probably) be done faster with binary searches (for large hulls)
        auto begin = std::find_if(outer_hull.begin(), outer_hull.end(), [&p](const indexed_segment& s) {
            return !is_left_of(p,s.first);
        });
        auto end = std::find_if(begin, outer_hull.end(), [&p](const indexed_segment& s) {
            return is_left_of(p,s.first);
        });
        for(auto it = begin; it != end; ++it) {
            triangles.push_back(Triangle{it->first.first, p, it->first.second,{it->second}});
            triangles[it->second].adjacent_indices.push_back(triangles.size()-1);
            if(it!=begin) {
                triangles[triangles.size()-1].adjacent_indices.push_back(triangles.size()-2);
                triangles[triangles.size()-2].adjacent_indices.push_back(triangles.size()-1);
            }
        }
        indexed_segment new_outer1(Segment(begin->first.first,p),end_tri), new_outer2(Segment(p,(end-1)->first.second),triangles.size()-1);
        if(begin==outer_hull.begin()) {
            auto begin2 = std::find_if(end, outer_hull.end(), [&p](const indexed_segment& s) {
                return !is_left_of(p,s.first);
            });
            
            if(begin2 != outer_hull.end()) {
                std::size_t end_tri2 = triangles.size();
                for(auto it = begin2; it != outer_hull.end(); ++it) {
                    triangles.push_back(Triangle{it->first.first, p, it->first.second,{it->second}});
                    triangles[it->second].adjacent_indices.push_back(triangles.size()-1);
                    if(it!=begin2) {
                        triangles[triangles.size()-1].adjacent_indices.push_back(triangles.size()-2);
                        triangles[triangles.size()-2].adjacent_indices.push_back(triangles.size()-1);
                    }
                }
                triangles.back().adjacent_indices.push_back(end_tri);
                triangles[end_tri].adjacent_indices.push_back(triangles.size()-1);
                new_outer1 = indexed_segment(Segment(begin2->first.first,p), end_tri2);
                outer_hull.erase(begin2, outer_hull.end());
            }
        }
        outer_hull.erase(begin, end);
        for(const auto& s : {new_outer1, new_outer2})
            outer_hull.insert(std::upper_bound(outer_hull.begin(), outer_hull.end(), s, outer_hull_pred),s);
    }

    //Step 9
    if(!skip_delaunay) {
        std::stack<std::pair<std::size_t,std::size_t>> L;
        for(std::size_t i=0; i < triangles.size(); ++i)
            for(const auto& a : triangles[i].adjacent_indices)
                if(i<a)
                    L.emplace(i,a);
        while(!L.empty()) {
            auto edge = L.top();
            L.pop();
            auto new_edges = flip_if_not_locally_delaunay<Point>(triangles, edge.first, edge.second);
            for(const auto& ne : new_edges)
                L.push(ne);
        }
    }
    std::vector<Segment> no_index_hull;
    no_index_hull.reserve(outer_hull.size());
    for(const auto& si : outer_hull)
        no_index_hull.push_back(si.first);
    return std::make_pair(triangles,no_index_hull);
}
