/*--------------------------------------------------------------------
DIVEMesh
Copyright 2008-2025 Hans Bihs

This file is part of DIVEMesh.

DIVEMesh is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Hans Bihs
--------------------------------------------------------------------*/

#include<iostream>
#include<math.h>
#include"iterators.h"
#include"looping.h"

#ifndef INCREMENT_H_
#define INCREMENT_H_

#define PI 3.14159265359

using namespace std;

class increment
{
	public:
	increment();
	virtual ~increment();
	static int i,j,k,aa,bb,cc,n,q;
	static int qn;
	static int xma,yma,zma;
    static int marge;
    static int margin;

	int conv(double);
	int equals(double,double);

};
#endif
