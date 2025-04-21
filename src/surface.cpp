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

#include"surface.h"
#include"lexer.h"
#include"dive.h"

surface::surface()
{
}

surface::~surface()
{
}

void surface::start(lexer* p, dive* a)
{
	mem_alloc(p,a);
    direction(p,a);
    makesurf(p,a);
    makesurf_plate(p,a);
    makesurfsolid(p,a);
}

void surface::mem_alloc(lexer* p, dive* a)
{
    cout<<"surface"<<endl;
    a->surfcount=0;

	surfnum=0;

    LOOP
    if(a->flag(i,j,k)==-20)
    a->flag(i,j,k)=-21;

	LOOP
    if(a->flag(i,j,k)>0)
    {
        if(a->flag(i-1,j,k)<0)
        ++surfnum;

        if(a->flag(i+1,j,k)<0)
        ++surfnum;

        if(a->flag(i,j-1,k)<0)
        ++surfnum;

		if(a->flag(i,j+1,k)<0)
        ++surfnum;

        if(a->flag(i,j,k-1)<0)
        ++surfnum;

        if(a->flag(i,j,k+1)<0)
        ++surfnum;
    }

    mem_alloc_plate(p,a);

	cout<<"surfnum: "<<surfnum<<"  surfnum_solid: "<<a->surfcount_solid<<endl;


	a->Iarray(a->surf,surfnum+a->surfcount_solid,5);
}

void surface::makesurf(lexer* p, dive* a)
{
//-------

    LOOP
    if(a->flag(i,j,k)>0)
    {
        if(a->flag(i-1,j,k)<0 && (p->C21<=1||i>0))
        {
        a->surf[a->surfcount][0]=i;
        a->surf[a->surfcount][1]=j;
        a->surf[a->surfcount][2]=k;
        a->surf[a->surfcount][3]=1;
        a->surf[a->surfcount][4]=fabs(a->flag(i-1,j,k));
        a->surfcount++;
        }

        if(a->flag(i+1,j,k)<0 && (p->C21<=1||i<p->knox-1))
        {
        a->surf[a->surfcount][0]=i;
        a->surf[a->surfcount][1]=j;
        a->surf[a->surfcount][2]=k;
        a->surf[a->surfcount][3]=4;
        a->surf[a->surfcount][4]=fabs(a->flag(i+1,j,k));
        a->surfcount++;
        }

        if(a->flag(i,j-1,k)<0 && (p->C22<=1||j>0))
        {
        a->surf[a->surfcount][0]=i;
        a->surf[a->surfcount][1]=j;
        a->surf[a->surfcount][2]=k;
        a->surf[a->surfcount][3]=3;
        a->surf[a->surfcount][4]=fabs(a->flag(i,j-1,k));
        a->surfcount++;
        }

		if(a->flag(i,j+1,k)<0 && (p->C22<=1||j<p->knoy-1))
        {
        a->surf[a->surfcount][0]=i;
        a->surf[a->surfcount][1]=j;
        a->surf[a->surfcount][2]=k;
        a->surf[a->surfcount][3]=2;
        a->surf[a->surfcount][4]=fabs(a->flag(i,j+1,k));
        a->surfcount++;
        }

        if(a->flag(i,j,k-1)<0 && (p->C23<=1||k>0))
        {
        a->surf[a->surfcount][0]=i;
        a->surf[a->surfcount][1]=j;
        a->surf[a->surfcount][2]=k;
        a->surf[a->surfcount][3]=5;
        a->surf[a->surfcount][4]=fabs(a->flag(i,j,k-1));
        a->surfcount++;
        }

        if(a->flag(i,j,k+1)<0 && (p->C23<=1||k<p->knoz-1))
        {
        a->surf[a->surfcount][0]=i;
        a->surf[a->surfcount][1]=j;
        a->surf[a->surfcount][2]=k;
        a->surf[a->surfcount][3]=6;
        a->surf[a->surfcount][4]=fabs(a->flag(i,j,k+1));
        a->surfcount++;
        }
    }

    cout<<"surfcount: "<<a->surfcount<<endl;


    for(i=0;i<a->surfcount;i++)
    {
    //cout<<" surface: "<<a->surf[i][0]<<" "<<a->surf[i][1]<<" "<<a->surf[i][2]<<" "<<a->surf[i][3]<<" "<<a->surf[i][4]<<endl;
    n=a->subgrid(a->surf[i][0],a->surf[i][1],a->surf[i][2]);
    a->wall[n]++;
    }

}
