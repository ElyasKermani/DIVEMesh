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
#include<fstream>
#include"lexer.h"
#include"increment.h"
#include"dive.h"

using namespace std;


#ifndef DECOMPOSITION_H_
#define DECOMPOSITION_H_

class decomp : public increment
{
public:
    decomp();
    virtual ~decomp();
    virtual void start(lexer*,dive*);
    
    void surfcount(lexer*,dive*);

private:
    void nodecalc(lexer*,dive*);
    void costfunc(lexer*,dive*);
	void partition_analyse(lexer*,dive*);
    void partition(lexer*,dive*);
    void partition_correction(lexer*, dive*);
	void partition_correct_x(lexer*, dive*);
	void partition_correct_y(lexer*, dive*);
	void partition_correct_z(lexer*, dive*);
	void partition_voidcheck(lexer*, dive*);
    void partition_manual(lexer*,dive*);
    void partition_planes(lexer*,dive*);
    void print_partition(lexer*, dive*);
    void rank(lexer*,dive*,int&,int&,int&,int&,int&,int&);

    void parasurface(lexer*,dive*);
    void paracosurface(lexer*,dive*);
    void cornersurface(lexer*,dive*);
    void paraslicesurface(lexer*,dive*);
    void paraslicecosurface(lexer*,dive*);
    void neighbors(lexer*,dive*);
    void cornercount(lexer*,dive*);
    void slicesurfcount(lexer*,dive*);
    void slicecornercount(lexer*,dive*);
    void knoxcalc(lexer*,dive*);
	void mem_alloc(lexer*,dive*);
    int partition_check(lexer*,dive*);
    
    void periodic_ini(lexer*,dive*);
    void periodic_surf(lexer*,dive*);
    void periodic_surfco(lexer*,dive*);
    void periodic_nb(lexer*,dive*);
    void periodic_count(lexer*,dive*);
    
    void decomp_vtp(lexer*,dive*);

    int numfactor;
    int active,nonactive,nodes,partqual,maxel;
    int *partcount;
    double factor, alpha;
    double cost, costmin;

    int posx,posy,posz;
    int count,pcount;
    int id;

    int jmin,jmax,**kmin,**kmax;
    int bcount;
	
	double xvar,yvar,zvar;
	double xvar_per,yvar_per,zvar_per;
	double xdiff,ydiff,zdiff;
	double xcross_m, ycross_m, zcross_m;
	int maindir;
	int *subcell;
	int *xcross,*ycross,*zcross;
	int *xcount,*ycount,*zcount;
	int xaverage,yaverage,zaverage;
	int xcount_sum,ycount_sum,zcount_sum,cellnum,paverage;
  
    
    char name[100],pname[100],epsvar[100];
    int n,iin,offset[100];
    float ffn;
    
    double **trix,**triy,**triz;
    int tricount;
	
	ofstream ddout;


};

#endif

