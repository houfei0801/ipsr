/*
Copyright (c) 2022, Fei Hou and Chiyu Wang, Institute of Software, Chinese Academy of Sciences.
All rights reserved.

The code can only be used for academic purpose, and cannot be used for commercial
purpose without written permission.

Redistribution and use in source for academic purpose, with or without modification,
are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list of
conditions and the following disclaimer.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO THE IMPLIED WARRANTIES
OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
TO, PROCUREMENT OF SUBSTITUTE  GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.
*/

#ifndef UTILITY_H
#define UTILITY_H

#include <vector>
#include <string>
#include <climits>
#include <fstream>
#include "PointStream.h"
#include "PointStreamData.h"

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
bool output_points_and_normals(const std::string &outFile, const std::vector<std::pair<Point<Real, Dim>, Normal<Real, Dim>>> &points_normals, const XForm<Real, Dim + 1> &iXForm)
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
		Point<Real, Dim> p = iXForm * points_normals[i].first;
		plyfile << p[0] << " " << p[1] << " " << p[2] << " ";
		plyfile << points_normals[i].second.normal[0] << " " << points_normals[i].second.normal[1] << " " << points_normals[i].second.normal[2] << std::endl;
	}
	plyfile.close();
	return true;
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
