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

#include"printer.h"
#include"increment.h"

class lexer;
class dive;

#ifndef PRINT_STL_H_
#define PRINT_STL_H_

using namespace std;

class print_stl :  public increment
{

public:
	print_stl(lexer*,dive*);
	virtual ~print_stl();
	virtual void solid_vtp(lexer*,dive*);
	virtual void solid_stl(lexer*,dive*);
    virtual void topo_vtp(lexer*,dive*);
	virtual void topo_stl(lexer*,dive*);

private:

    char name[100],pname[100],epsvar[100];
    int n,iin,offset[100];
    float ffn;
    double ddn;
};

#endif


