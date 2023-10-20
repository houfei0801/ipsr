/*
Copyright (c) 2022, Fei Hou and Chiyu Wang, Institute of Software, Chinese Academy of Sciences.
All rights reserved.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#ifndef UTILITY_H
#define UTILITY_H

#include <vector>
#include <string>
#include <climits>
#include <fstream>
#include "kdtree.h"
#include "PointStream.h"
#include "PointStreamData.h"

template <class Real, unsigned int Dim>
void transform(std::vector<std::pair<Point<Real, Dim>, Normal<Real, Dim>>> &points_normals, const XForm<Real, Dim + 1> &iXForm)
{
	for (size_t i = 0; i < points_normals.size(); ++i)
	{
		points_normals[i].first = iXForm * points_normals[i].first;
	}
}

template <class Real, unsigned int Dim>
void ply_reader(const std::string &file, std::vector<std::pair<Point<Real, Dim>, Normal<Real, Dim>>> &points_normals)
{
	PLYInputPointStream<Real, Dim> ply(file.c_str());
	Normal<Real, Dim> n(Point<Real, Dim>(1, 0, 0));
	Point<Real, Dim> p;
	while (ply.nextPoint(p))
		points_normals.push_back(std::make_pair(p, n));
}

template <class Real, unsigned int Dim>
bool output_ply(const std::string &outFile, const std::pair<std::vector<Point<Real, Dim>>, std::vector<std::vector<int>>> &mesh, const XForm<Real, Dim + 1> &iXForm)
{
	const std::vector<Point<Real, Dim>> &points = mesh.first;
	const std::vector<std::vector<int>> &faces = mesh.second;

	std::ofstream plyfile;
	plyfile.open(outFile, std::ofstream::out);
	if (!plyfile)
	{
		printf("Cannot save result file %s\n", outFile.c_str());
		return false;
	}
	printf("writing to %s\n", outFile.c_str());

	plyfile << "ply\nformat ascii 1.0\n";
	plyfile << "element vertex " << points.size() << std::endl;
	plyfile << "property float x" << std::endl
		<< "property float y" << std::endl
		<< "property float z" << std::endl;
	plyfile << "element face " << faces.size() << std::endl;
	plyfile << "property list uchar int vertex_index" << std::endl;
	plyfile << "end_header" << std::endl;

	for (size_t i = 0; i < points.size(); ++i)
	{
		Point<Real, Dim> p = iXForm * points[i];
		plyfile << p[0] << " " << p[1] << " " << p[2] << std::endl;
	}

	for (size_t i = 0; i < faces.size(); ++i)
	{
		plyfile << faces[i].size();
		for (size_t j = 0; j < faces[i].size(); ++j)
			plyfile << " " << faces[i][j];
		plyfile << std::endl;
	}
	plyfile.close();
	return true;
}

template <class Real, unsigned int Dim>
bool output_sample_points_and_normals(const std::string &outFile, const std::vector<std::pair<Point<Real, Dim>, Normal<Real, Dim>>> &points_normals, const XForm<Real, Dim + 1> &iXForm)
{
	std::ofstream plyfile;
	plyfile.open(outFile, std::ofstream::out);
	if (!plyfile)
	{
		printf("Cannot save result file %s\n", outFile.c_str());
		return false;
	}
	printf("writing to %s\n", outFile.c_str());

	plyfile << "ply\nformat ascii 1.0\n";
	plyfile << "element vertex " << points_normals.size() << std::endl;
	plyfile << "property float x" << std::endl
		<< "property float y" << std::endl
		<< "property float z" << std::endl;
	plyfile << "property float nx" << std::endl
		<< "property float ny" << std::endl
		<< "property float nz" << std::endl;
	plyfile << "element face 0" << std::endl;
	plyfile << "property list uchar int vertex_index" << std::endl;
	plyfile << "end_header" << std::endl;

	for (size_t i = 0; i < points_normals.size(); ++i)
	{
		auto p = iXForm * points_normals[i].first;
		plyfile << p[0] << " " << p[1] << " " << p[2] << " ";
		plyfile << points_normals[i].second.normal[0] << " " << points_normals[i].second.normal[1] << " " << points_normals[i].second.normal[2] << std::endl;
	}
	plyfile.close();
	return true;
}

template <class Real, unsigned int Dim>
bool output_all_points_and_normals(const std::string &outFile, const std::string &input_name, const std::vector<std::pair<Point<Real, Dim>, Normal<Real, Dim>>> &points_normals, const kdt::KDTree<kdt::KDTreePoint> &tree, const XForm<Real, Dim + 1> &iXForm)
{
	std::vector<std::pair<Point<Real, Dim>, Normal<Real, Dim>>> points_normals_all;
	ply_reader<Real, Dim>(input_name, points_normals_all);
	auto inv_iXForm = iXForm.inverse();
	for (size_t i = 0; i < points_normals_all.size(); ++i)
	{
		auto c = inv_iXForm * points_normals_all[i].first;
		std::array<Real, Dim> a{ c[0], c[1], c[2] };
		int n = tree.nnSearch(kdt::KDTreePoint(a));
		points_normals_all[i].second = points_normals[n].second;
	}

	std::ofstream plyfile;
	plyfile.open(outFile, std::ofstream::out);
	if (!plyfile)
	{
		printf("Cannot save result file %s\n", outFile.c_str());
		return false;
	}
	printf("writing to %s\n", outFile.c_str());

	plyfile << "ply\nformat ascii 1.0\n";
	plyfile << "element vertex " << points_normals_all.size() << std::endl;
	plyfile << "property float x" << std::endl
		<< "property float y" << std::endl
		<< "property float z" << std::endl;
	plyfile << "property float nx" << std::endl
		<< "property float ny" << std::endl
		<< "property float nz" << std::endl;
	plyfile << "element face 0" << std::endl;
	plyfile << "property list uchar int vertex_index" << std::endl;
	plyfile << "end_header" << std::endl;

	for (size_t i = 0; i < points_normals_all.size(); ++i)
	{
		const auto& p =  points_normals_all[i].first;
		plyfile << p[0] << " " << p[1] << " " << p[2] << " ";
		plyfile << points_normals_all[i].second.normal[0] << " " << points_normals_all[i].second.normal[1] << " " << points_normals_all[i].second.normal[2] << std::endl;
	}
	plyfile.close();
	return true;
}

template <class Real, int Dim>
bool operator==(const Normal<Real, Dim> &n1, const Normal<Real, Dim> &n2)
{
	for (int i = 0; i < Dim; ++i)
		if (n1.normal[i] != n2.normal[i])
			return false;
	return true;
}

template <class Real, unsigned int Dim>
void normalize(Normal<Real, Dim> &n)
{
	Real len = 0;
	for (unsigned int i = 0; i < Dim; ++i)
		len += n.normal[i] * n.normal[i];
	if (len != 0)
	{
		len = sqrt(len);
		for (unsigned int i = 0; i < Dim; ++i)
			n.normal[i] /= len;
	}
}

inline std::vector<std::string> split(const std::string &s, char c = ' ')
{
	std::vector<std::string> str;
	unsigned pos = 0;
	while (pos < s.size())
	{
		while (pos < s.size() && s[pos] == c)
			++pos;

		unsigned int end = pos;

		do
		{
			++end;
		} while (end < s.size() && s[end] != c);

		if (pos < s.size())
			str.push_back(s.substr(pos, end - pos));
		pos = end;
	}

	return str;
}

inline bool valid_parameter(long v)
{
	return v > 0 && v < INT_MAX;
}

#endif
