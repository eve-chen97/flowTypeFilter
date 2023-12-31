#include <iostream>
#include <list>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>


int main()
{
	typedef boost::geometry::model::d2::point_xy<double> point_type;
	typedef boost::geometry::model::polygon<point_type> polygon_type;

	polygon_type poly;
	boost::geometry::read_wkt(
		"POLYGON((2 1.3,2.4 1.7,2.8 1.8,3.4 1.2,3.7 1.6,3.4 2,4.1 3,5.3 2.6,5.4 1.2,4.9 0.8,2.9 0.7,2 1.3)"
		"(4.0 2.0, 4.2 1.4, 4.8 1.9, 4.4 2.2, 4.0 2.0))", poly);

	point_type p(4, 1);

	std::cout << "within: " << (boost::geometry::within(p, poly) ? "yes" : "no") << std::endl;

	return 0;
}