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
#include <fstream>

void lexer::read_control_preread()
{
	char c;
	int numint,numval;
    int count;
	numlinsurf=0;

	cout<<"read control"<<endl;

	ifstream control("control.txt", ios_base::in);

	if(!control)
	{
		cout<<"no 'control.txt' file"<<endl<<endl;
		exit(1);
	}
    
    count=0;
	while(!control.eof())
	{
	    control>>c;

	if (c == '/') 
	{
		control.ignore(1000, '\n');
	}
	else
	{	
		switch(c)
		{
        case 'B': control>>numint;
				switch(numint)
				{
                case 3: control>>B3_dx>>B3_dy;
                          B3=1;
                          clear(c,numint);
						 break;
                case 4: control>>B4;
                          clear(c,numint);
						 break;
                case 5: control>>B5;
                          clear(c,numint);
						 break;
                case 10: ++B10;
                          clear(c,numint);
						 break;
                case 22: ++B22;
						 clear(c,numint);
						 break;
                case 31: ++B31;
						 clear(c,numint);
						 break;
                case 32: ++B32;
						 clear(c,numint);
						 break;
                case 33: ++B33;
						 clear(c,numint);
						 break;
                case 101: control>>B101;
						 clear(c,numint);
						 break;
                case 102: control>>B102;
						 clear(c,numint);
						 break;
                case 103: control>>B103;
						 clear(c,numint);
						 break;
                case 111: control>>B111;
						 clear(c,numint);
						 break;
                case 112: control>>B112;
						 clear(c,numint);
						 break;
                case 113: control>>B113;
						 clear(c,numint);
						 break;
                case 114: control>>B114_x;
                        B114=1;
						 clear(c,numint);
						 break;
                case 115: control>>B115_y;
                        B115=1;
						 clear(c,numint);
						 break;
                case 116: control>>B116_z;
                        B116=1;
						 clear(c,numint);
						 break;
                case 121: control>>B121_N1>>B121_x1>>B121_N2>>B121_x2>>B121_N3;
                        B121=1;
						 clear(c,numint);
                        break;
                case 122: control>>B122_N1>>B122_y1>>B122_N2>>B122_y2>>B122_N3;
                        B122=1;
						 clear(c,numint);
                        break;
                case 123: control>>B123_N1>>B123_z1>>B123_N2>>B123_z2>>B123_N3;
                        B123=1;
						 clear(c,numint);
                        break;
                case 124: control>>B124_N1>>B124_x1>>B124_f1>>B124_N2>>B124_x2>>B124_f2>>B124_N3>>B124_x3>>B124_f3;
                        B124=1;
						 clear(c,numint);
						 break;
                case 125: control>>B125_N1>>B125_y1>>B125_f1>>B125_N2>>B125_y2>>B125_f2>>B125_N3>>B125_y3>>B125_f3;
                        B125=1;
						 clear(c,numint);
						 break;
                case 126: control>>B126_N1>>B126_z1>>B126_f1>>B126_N2>>B126_z2>>B126_f2>>B126_N3>>B126_z3>>B126_f3;
                        B126=1;
						 clear(c,numint);
						 break;
                case 127: control>>B127_dx_min>>B127_dx_max>>B127_pf>>B127_df>>B127_r;
                        B127=1;
						 clear(c,numint);
						 break;
                case 128: control>>B128_dx_min>>B128_dx_max>>B128_pf>>B128_df>>B128_r;
                        B128=1;
						 clear(c,numint);
						 break;   
                case 129: control>>B129_dx_min>>B129_dx_max>>B129_pf>>B129_df>>B129_r;
                        B129=1;
						 clear(c,numint);
						 break;  
				case 130: control>>B130;
						 clear(c,numint);
						 break;                      
				}
				break;
                
       
				
        case 'S': 
        case 'O': 
        control>>numint;
				switch(numint)
				{
				case 1: control>>S1;
						 clear(c,numint);
						 break;
				case 2: control>>S2;
						 clear(c,numint);
						 break;
				case 3: control>>S3_xs>>S3_xe>>S3_ys>>S3_ye>>S3_zs>>S3_ze;
						S3=1;
						 clear(c,numint);
						 break;
				case 4: control>>S4;
						 clear(c,numint);
						 break;
				case 5: control>>S5_x>>S5_y>>S5_z>>S5_phi>>S5_theta>>S5_psi;
                          S5=1;
						 clear(c,numint);
						 break;
				case 7: control>>S7_dx>>S7_dy>>S7_dz;
						S7=1;
						 clear(c,numint);
						 break;
				case 8: control>>S8;
						 clear(c,numint);
						 break;
				
				}
				break;
                
        case 'T': 
        control>>numint;
				switch(numint)
				{
				case 1: control>>T1;
						 clear(c,numint);
						 break;
				case 2: control>>T2;
						 clear(c,numint);
						 break;
				case 3: control>>T3_xs>>T3_xe>>T3_ys>>T3_ye>>T3_zs>>T3_ze;
						T3=1;
						 clear(c,numint);
						 break;
				case 4: control>>T4;
						 clear(c,numint);
						 break;
				case 5: control>>T5_x>>T5_y>>T5_z>>T5_phi>>T5_theta>>T5_psi;
                          T5=1;
						 clear(c,numint);
						 break;
				case 7: control>>T7_dx>>T7_dy>>T7_dz;
						T7=1;
						 clear(c,numint);
						 break;
				case 8: control>>T8;
						 clear(c,numint);
						 break;
				}
				break;


		}
        ++count;
	}
        if(count>1e7)
        {
        cout<<endl;
        cout<<"!!! missing input parameter in control.txt !!!"<<endl<<endl;
        cout<<"!!! please check the DIVEMesh User Guide !!!"<<endl<<endl<<endl<<endl;
        
        exit(1);
        }
	}
	control.close();
	control.clear();
	
}