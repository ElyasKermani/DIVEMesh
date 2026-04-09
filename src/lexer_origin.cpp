/*--------------------------------------------------------------------
DIVEMesh
Copyright 2008-2026 Hans Bihs

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

#include"lexer.h"

void lexer::origin_calc()
{
    // xmin,xmax
    xmin=+1.0e19;
    ymin=+1.0e19;
    zmin=+1.0e19;

    xmax=-1.0e19;
    ymax=-1.0e19;
    zmax=-1.0e19;

    for(qn=0;qn<B10;++qn)
    {
    xmin=MIN(xmin,B10_xs[qn]);
    ymin=MIN(ymin,B10_ys[qn]);
    zmin=MIN(zmin,B10_zs[qn]);
    xmin=MIN(xmin,B10_xe[qn]);
    ymin=MIN(ymin,B10_ye[qn]);
    zmin=MIN(zmin,B10_ze[qn]);

    xmax=MAX(xmax,B10_xs[qn]);
    ymax=MAX(ymax,B10_ys[qn]);
    zmax=MAX(zmax,B10_zs[qn]);
    xmax=MAX(xmax,B10_xe[qn]);
    ymax=MAX(ymax,B10_ye[qn]);
    zmax=MAX(zmax,B10_ze[qn]);
    }

    for(qn=0;qn<B22;++qn)
    {
    xmin=MIN(xmin,B22_xm[qn]-B22_r[qn])-dx;
    ymin=MIN(ymin,B22_ym[qn]-B22_r[qn])-dx;
    zmin=MIN(zmin,B22_zm[qn]-B22_r[qn])-dx;

    xmax=MAX(xmax,B22_xm[qn]+B22_r[qn])+dx;
    ymax=MAX(ymax,B22_ym[qn]+B22_r[qn])+dx;
    zmax=MAX(zmax,B22_zm[qn]+B22_r[qn])+dx;
    }

    for(qn=0;qn<B31;++qn)
    {
    xmin=MIN(xmin,B31_xs[qn]);
    xmin=MIN(xmin,B31_xe[qn]);
    ymin=MIN(ymin,B31_ym[qn]-B31_r[qn])-dx;
    zmin=MIN(zmin,B31_zm[qn]-B31_r[qn])-dx;

    xmax=MAX(xmax,B31_xs[qn]);
    xmax=MAX(xmax,B31_xe[qn]);
    ymax=MAX(ymax,B31_ym[qn]+B31_r[qn])+dx;
    zmax=MAX(zmax,B31_zm[qn]+B31_r[qn])+dx;
    }

    for(qn=0;qn<B32;++qn)
    {
    xmin=MIN(xmin,B32_xm[qn]-B32_r[qn])-dx;
    ymin=MIN(ymin,B32_ys[qn]);
    ymin=MIN(ymin,B32_ye[qn]);
    zmin=MIN(zmin,B32_zm[qn]-B32_r[qn])-dx;

    xmax=MAX(xmax,B32_xm[qn]+B32_r[qn])+dx;
    ymax=MAX(ymax,B32_ys[qn]);
    ymax=MAX(ymax,B32_ye[qn]);
    zmax=MAX(zmax,B32_zm[qn]+B32_r[qn])+dx;
    }

    for(qn=0;qn<B33;++qn)
    {
    xmin=MIN(xmin,B33_xm[qn]-B33_r[qn])-dx;
    ymin=MIN(ymin,B33_ym[qn]-B33_r[qn])-dx;
    zmin=MIN(zmin,B33_zs[qn]);
    zmin=MIN(zmin,B33_ze[qn]);

    xmax=MAX(xmax,B33_xm[qn]+B33_r[qn])+dx;
    ymax=MAX(ymax,B33_ym[qn]+B33_r[qn])+dx;
    zmax=MAX(zmax,B33_zs[qn]);
    zmax=MAX(zmax,B33_ze[qn]);
    }
    
    if(S1==1 && S2==1)
    {
    xmin=xs_stl;
    ymin=ys_stl;
    zmin=zs_stl;
    xmax=xe_stl;
    ymax=ye_stl;
    zmax=ze_stl;
	}
    
    // origin
    global_orig_x = xmin + B3_dx;
    global_orig_y = ymin + B3_dy;
    
    
}