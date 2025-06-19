// $Source$
//--------------------------------------------------------------------------------
// global
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/04/25
//
/** @file global.cpp
 *  @brief Implements functions for loading Earth orientation parameters, gravity models, planetary ephemeris coefficients, and observation data.
 *
 *  @author Jarein López Ruiz
 *  @bug No known bugs.
 */
//--------------------------------------------------------------------------------

#include "..\include\global.hpp"

Param AuxParam;
Matrix eopdata;
Matrix Cnm;
Matrix Snm;
int n_eqn;

void eop19620101(int c){
	eopdata=zeros(13,c);
	
	FILE *fp = fopen("../data/eop19620101.txt","r");
	if(fp == NULL){
		printf("Fail open eop19620101.txt file\n");
		exit(EXIT_FAILURE);
	}
	
	for(int j=1;j<=c;j++){
		fscanf(fp,"%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",&(eopdata(1,j)),
					&(eopdata(2,j)),&(eopdata(3,j)),&(eopdata(4,j)),&(eopdata(5,j)),
					&(eopdata(6,j)),&(eopdata(7,j)),&(eopdata(8,j)),&(eopdata(9,j)),
					&(eopdata(10,j)),&(eopdata(11,j)),&(eopdata(12,j)),&(eopdata(13,j)));
	}
	
	fclose(fp);
}

void GGM03S(int n){
	Cnm=zeros(n,n);
	Snm=zeros(n,n);
	
	FILE *fp = fopen("../data/GGM03S.txt","r");
	if(fp == NULL){
		printf("Fail open GGM03S.txt file\n");
		exit(EXIT_FAILURE);
	}
	
	double aux;
	for(int i=1;i<=n;i++){
		for (int j=1;j<=i;j++){
			fscanf(fp,"%lf%lf%lf%lf%lf%lf",&aux,&aux,&(Cnm(i,j)),&(Snm(i,j)),&aux,&aux);
		}
	}
	
	fclose(fp);
}

Matrix PC;

void DE430Coeff(int row, int col) {
    PC = zeros(row, col);

    FILE *fp = fopen("../data/DE430Coeff.txt", "r");
    if (fp == NULL) {
        printf("Fail open DE430Coeff.txt file\n");
        exit(EXIT_FAILURE);
    }

    double aux;
    int total_columns = 1020;

    for (int i = 1; i <= row; i++) {
        for (int j = 1; j <= col; j++) {
            fscanf(fp, "%lf", &PC(i, j));
        }
        for (int k = col + 1; k <= total_columns; k++) {
            fscanf(fp, "%lf", &aux);
        }
    }

    fclose(fp);
}

Matrix obs;

void GEOS3(int f){
	obs=zeros(f,4);
	
	FILE *fid = fopen("../data/GEOS3.txt","r");
	if(fid==NULL){
		printf("Fail open GEOS3.txt file\n");
		exit(EXIT_FAILURE);
	}
	
	int Y,MO,D,H,MI,S;
	double AZ,EL,DIST;
	char line[55],y[5],mo[3],d[3],h[3],mi[3],s[3],az[9],el[8],dist[11];
	for(int i = 1; i <= f; i++) {
    char line[128];
    if (fgets(line, sizeof(line), fid) == NULL) break;

    int Y, MO, D, H, MI;
    double S, AZ, EL, DIST;

    // Formato adaptado al fichero real (espacios y valores con decimales)
    sscanf(line, "%d/%d/%d %d:%d:%lf %lf %lf %lf",
           &Y, &MO, &D, &H, &MI, &S, &AZ, &EL, &DIST);

    obs(i,1) = Mjday(Y, MO, D, H, MI, S);
    obs(i,2) = SAT_Const::Rad * AZ;
    obs(i,3) = SAT_Const::Rad * EL;
    obs(i,4) = 1e3 * DIST;  // pasa de km a m
}


	fclose(fid);
}
