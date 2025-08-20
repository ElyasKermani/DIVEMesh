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

#include"geodat.h"
#include"dive.h"
#include"lexer.h"
#include"interpolation.h"

void geodat::geo_patch(lexer *p, dive *a, field2d& bed)
{
    p->Darray(gpf,kx+2*dd+1,ky+2*dd+1);
    
    // read geo_patch
    geo_patch_read(p,a);
    
    pipol->start(p,a,numGP,GP_x,GP_y,GP_z,XC,YC,kx,ky,gpf);
    
    // correct topof
    double xc,yc,val;

    cout<<"patch bed "<<endl;
    
    XYLOOP
    {
    xc = p->XP[IP];
    yc = p->YP[JP];
    
    val = ccipol(p,gpf,xc,yc); 

    //cout<<" PROLONG: "<<topof[i+3][j+3]<<" val: "<<val<<endl;  

    bed(i,j) = MAX(val,bed(i,j));    
    }
    
}



void geodat::geo_patch_read(lexer *p, dive *a)
{
// read geo_patch

        cout<<"open geo_patch.dat and count entries"<<endl;
        double val;
        char cval;
        ifstream geo("geo_patch.dat", ios_base::in);

        while(!geo.eof())
        {
        if(p->G19==0)
        geo>>val>>val>>val;

        if(p->G19==1)
        geo>>cval>>val>>val>>val;

        ++countGP;
        }
        
        numGP=countGP-1;
        cout<<"> geo_patch entries: "<<numGP<<endl;
        geo.close();

        GP_x = new double[countGP];
        GP_y = new double[countGP];
        GP_z = new double[countGP];

        
        cout<<"read geo_patch.dat"<<endl;
        
        geo.open("geo_patch.dat", ios_base::in);

        countGP=0;
        while(!geo.eof()&&countGP<numGP)
        {
        if(p->G19==0)
        geo>>GP_x[countGP]>>GP_y[countGP]>>GP_z[countGP];

        if(p->G19==1)
        geo>>cval>>GP_x[countGP]>>GP_y[countGP]>>GP_z[countGP];
        ++countGP;
        }

        geo.close();
		
		if(p->G13>0)
		{   

			
			double xval,yval;
			for(n=0;n<numGP;++n)
			{
			xval = p->G14_x + (GP_x[n]-p->G14_x)*cos(p->G13_phi) - (GP_y[n]-p->G14_y)*sin(p->G13_phi);
			yval = p->G14_y + (GP_x[n]-p->G14_x)*sin(p->G13_phi) + (GP_y[n]-p->G14_y)*cos(p->G13_phi);

			GP_x[n] = xval;
			GP_y[n] = yval;
			}
		}
        

}