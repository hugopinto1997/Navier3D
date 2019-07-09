void showMatrix(Matrix K){
    for(int i=0;i<K.at(0).size();i++){
        cout << "[\t";
        for(int j=0;j<K.size();j++){
            cout << K.at(i).at(j) << "\t";
        }
        cout << "]\n";
    }
}

void showKs(vector<Matrix> Ks){
    for(int i=0;i<Ks.size();i++){
        cout << "K del elemento "<< i+1 << ":\n";
        showMatrix(Ks.at(i));
        cout << "*************************************\n";
    }
}

void showVector(Vector b){
    cout << "[\t";
    for(int i=0;i<b.size();i++){
        cout << b.at(i) << "\t";
    }
    cout << "]\n";
}

void showbs(vector<Vector> bs){
    for(int i=0;i<bs.size();i++){
        cout << "b del elemento "<< i+1 << ":\n";
        showVector(bs.at(i));
        cout << "*************************************\n";
    }
}

node selectNode(int i, element e,mesh &m){
	node n;
	switch(i){
		case 1: n = m.getNode(e.getNode1()-1); break;
		case 2: n = m.getNode(e.getNode2()-1); break;
		case 3: n = m.getNode(e.getNode3()-1); break;
        case 4: n = m.getNode(e.getNode4()-1); break;
	}
	return n;
}

float selectCoord(int c, node n){
	float v;
	switch(c){
		case EQUIS: v = n.getX(); break;
		case YE: v = n.getY(); break;
        case ZETA: v = n.getZ(); break;
	}
	return v;
}

float calcularTenedor(element e, int coord, int i, int j,mesh &m){
	node n1=selectNode(i,e,m),n2=selectNode(j,e,m);

	return selectCoord(coord,n1) - selectCoord(coord,n2);
}

float calculateLocalD(int i,mesh m){
    Matrix matriz;
    Vector row1, row2, row3;

    element e = m.getElement(i);
    node n1 = m.getNode(e.getNode1()-1);
    node n2 = m.getNode(e.getNode2()-1);
    node n3 = m.getNode(e.getNode3()-1);

	row1.push_back(calcularTenedor(e,EQUIS,2,1,m)); 
    row1.push_back(calcularTenedor(e,YE,2,1,m)); 
    row1.push_back(calcularTenedor(e,ZETA,2,1,m));

	row2.push_back(calcularTenedor(e,EQUIS,3,1,m)); 
    row2.push_back(calcularTenedor(e,YE,3,1,m));
    row2.push_back(calcularTenedor(e,ZETA,3,1,m));

    row3.push_back(calcularTenedor(e,EQUIS,4,1,m)); 
    row3.push_back(calcularTenedor(e,YE,4,1,m));
    row3.push_back(calcularTenedor(e,ZETA,4,1,m));

	matriz.push_back(row1); matriz.push_back(row2); matriz.push_back(row3);

    return determinant(matriz);
}

float calculateMagnitude(float v1, float v2, float v3){
    return sqrt(pow(v1,2)+pow(v2,2)+pow(v3, 2));
}

float calculateLocalVolumen(int i,mesh m){

    double Ve, u, v, w, U, V, W, a, b, c, d, X, x, Y, y, Z, z;
    
    element e = m.getElement(i);
    node n1 = m.getNode(e.getNode1()-1);
    node n2 = m.getNode(e.getNode2()-1);
    node n3 = m.getNode(e.getNode3()-1);
    node n4 = m.getNode(e.getNode4()-1);

    U = calculateMagnitude(n2.getX()-n1.getX(), n2.getY()-n1.getY(), n2.getZ()-n1.getZ() );
    V = calculateMagnitude(n3.getX()-n2.getX(), n3.getY()-n2.getY(), n3.getZ()-n2.getZ() );
    W = calculateMagnitude(n3.getX()-n1.getX(), n3.getY()-n1.getY(), n3.getZ()-n1.getZ() );

    u = calculateMagnitude(n4.getX()-n3.getX(), n4.getY()-n3.getY(), n4.getZ()-n3.getZ() );
    v = calculateMagnitude(n4.getX()-n1.getX(), n4.getY()-n1.getY(), n4.getZ()-n1.getZ() );
    w = calculateMagnitude(n4.getX()-n2.getX(), n4.getY()-n2.getY(), n4.getZ()-n2.getZ() );

    X = (w-U+v)*(U+v+w);
    x = (U-v+w)*(v-w+U);

    Y = (u-V+w)*(V+w+u);
    y = (V-w+u)*(w-u+V);

    Z = (v-W+u)*(W+u+v);
    z = (W-u+v)*(u-v+W);

    a = sqrt(x*Y*Z);
    b = sqrt(y*Z*X);
    c = sqrt(z*X*Y);
    d = sqrt(x*y*z);

    Ve = sqrt( (-a+b+c+d)*(a-b+c+d)*(a+b-c+d)*(a+b+c-d) ) / (192*u*v*w);
    
    return Ve;

}

float OperarRestaTenedor(float a, float b, float c, float d){
    return (a*b)-(c*d);
}

void calculateLocalA(int i,Matrix &A,mesh m){
    zeroes(A,3);
    element e = m.getElement(i);
    node n1 = m.getNode(e.getNode1()-1);
    node n2 = m.getNode(e.getNode2()-1);
    node n3 = m.getNode(e.getNode3()-1);

    A.at(0).at(0) = OperarRestaTenedor(calcularTenedor(e,YE,3,1,m),calcularTenedor(e,ZETA,4,1,m),calcularTenedor(e,YE,4,1,m),calcularTenedor(e,ZETA,3,1,m));
    A.at(0).at(1) = OperarRestaTenedor(calcularTenedor(e,YE,4,1,m),calcularTenedor(e,ZETA,2,1,m),calcularTenedor(e,YE,2,1,m),calcularTenedor(e,ZETA,4,1,m));
    A.at(0).at(2) = OperarRestaTenedor(calcularTenedor(e,YE,2,1,m),calcularTenedor(e,ZETA,3,1,m),calcularTenedor(e,YE,3,1,m),calcularTenedor(e,ZETA,2,1,m));

    A.at(1).at(0) = OperarRestaTenedor(calcularTenedor(e,EQUIS,4,1,m),calcularTenedor(e,ZETA,3,1,m),calcularTenedor(e,EQUIS,3,1,m),calcularTenedor(e,ZETA,4,1,m));
    A.at(1).at(1) = OperarRestaTenedor(calcularTenedor(e,EQUIS,2,1,m),calcularTenedor(e,ZETA,4,1,m),calcularTenedor(e,EQUIS,4,1,m),calcularTenedor(e,ZETA,2,1,m));
    A.at(1).at(2) = OperarRestaTenedor(calcularTenedor(e,EQUIS,3,1,m),calcularTenedor(e,ZETA,2,1,m),calcularTenedor(e,EQUIS,2,1,m),calcularTenedor(e,ZETA,3,1,m));

    A.at(2).at(0) = OperarRestaTenedor(calcularTenedor(e,EQUIS,3,1,m),calcularTenedor(e,YE,4,1,m),calcularTenedor(e,EQUIS,4,1,m),calcularTenedor(e,YE,3,1,m));
    A.at(2).at(1) = OperarRestaTenedor(calcularTenedor(e,EQUIS,4,1,m),calcularTenedor(e,YE,2,1,m),calcularTenedor(e,EQUIS,2,1,m),calcularTenedor(e,YE,4,1,m));
    A.at(2).at(2) = OperarRestaTenedor(calcularTenedor(e,EQUIS,2,1,m),calcularTenedor(e,YE,3,1,m),calcularTenedor(e,EQUIS,3,1,m),calcularTenedor(e,YE,2,1,m));

}

void calculateB(Matrix &B){
    zeroes(B,3,12);

    B.at(0).at(0) = -1;
    B.at(0).at(1) = 1; 
    B.at(0).at(2) = 0; 
    B.at(0).at(3) = 0; 
    B.at(0).at(4) = -1; 
    B.at(0).at(5) = 1;
    B.at(0).at(6) = 0; 
    B.at(0).at(7) = 0; 
    B.at(0).at(8) = -1; 
    B.at(0).at(9) = 1; 
    B.at(0).at(10) = 0; 
    B.at(0).at(11) = 0;

    B.at(1).at(0) = -1;
    B.at(1).at(1) = 0; 
    B.at(1).at(2) = 1; 
    B.at(1).at(3) = 0; 
    B.at(1).at(4) = -1; 
    B.at(1).at(5) = 0;
    B.at(1).at(6) = 1; 
    B.at(1).at(7) = 0; 
    B.at(1).at(8) = -1; 
    B.at(1).at(9) = 0; 
    B.at(1).at(10) = 1; 
    B.at(1).at(11) = 0;

    B.at(2).at(0) = -1;
    B.at(2).at(1) = 0; 
    B.at(2).at(2) = 0; 
    B.at(2).at(3) = 1; 
    B.at(2).at(4) = -1; 
    B.at(2).at(5) = 0;
    B.at(2).at(6) = 0; 
    B.at(2).at(7) = 1; 
    B.at(2).at(8) = -1; 
    B.at(2).at(9) = 0; 
    B.at(2).at(10) = 0; 
    B.at(2).at(11) = 1;
}

void calculateC(Matrix &C){
    zeroes(C,3,4);
    C.at(0).at(0) = -1; C.at(0).at(1) = 1; C.at(0).at(2) = 0; C.at(0).at(3) = 0;
    C.at(1).at(0) = -1; C.at(1).at(1) = 0; C.at(1).at(2) = 1; C.at(1).at(3) = 0;
    C.at(2).at(0) = -1; C.at(2).at(1) = 0; C.at(2).at(2) = 0; C.at(2).at(3) = 1;
}

void ubicarSubMatriz(Matrix &K,int fi,int ff,int ci,int cf,Matrix M){
    int n = 0, m= 0;
    for(int i=fi;i<=ff;i++){
        for(int j=ci;j<=cf;j++){
            K.at(i).at(j) = M.at(n).at(m);
            m++;
        }
        n++; m = 0;
    }
}

void calculateEStar(Matrix& m){
	zeroes(m,12,3);

	m.at(0).at(0) = 0;          m.at(0).at(1) = 0;          m.at(0).at(2) = 0;
    m.at(1).at(0) = 0.5;        m.at(1).at(1) = 0;          m.at(1).at(2) = 0;
    m.at(2).at(0) = 0.5;        m.at(2).at(1) = 0;          m.at(2).at(2) = 0;
    m.at(3).at(0) = 0.5;        m.at(3).at(1) = 0;          m.at(3).at(2) = 0;
    m.at(4).at(0) = 0;          m.at(4).at(1) = 0;          m.at(4).at(2) = 0;
    m.at(5).at(0) = 0;          m.at(5).at(1) = 0.5;        m.at(5).at(2) = 0;
    m.at(6).at(0) = 0;          m.at(6).at(1) = 0.5;        m.at(6).at(2) = 0;
    m.at(7).at(0) = 0;          m.at(7).at(1) = 0.5;        m.at(7).at(2) = 0;
    m.at(8).at(0) = 0;          m.at(8).at(1) = 0;          m.at(8).at(2) = 0;
    m.at(9).at(0) = 0;          m.at(9).at(1) = 0;          m.at(9).at(2) = 0.5;
    m.at(10).at(0) = 0;         m.at(10).at(1) = 0;         m.at(10).at(2) = 0.5;
    m.at(11).at(0) = 0;         m.at(11).at(1) = 0;         m.at(11).at(2) = 0.5;
}

float calculateLocalJ(int i,mesh m){
    Matrix matriz;
    Vector row1, row2, row3;

    element e = m.getElement(i);
    node n1 = m.getNode(e.getNode1()-1);
    node n2 = m.getNode(e.getNode2()-1);
    node n3 = m.getNode(e.getNode3()-1);

	row1.push_back(calcularTenedor(e,EQUIS,2,1,m)); 
    row1.push_back(calcularTenedor(e,EQUIS,3,1,m));
    row1.push_back(calcularTenedor(e,EQUIS,4,1,m));
	
    row2.push_back(calcularTenedor(e,YE,2,1,m)); 
    row2.push_back(calcularTenedor(e,YE,3,1,m));
    row2.push_back(calcularTenedor(e,YE,4,1,m));

    row3.push_back(calcularTenedor(e,ZETA,2,1,m)); 
    row3.push_back(calcularTenedor(e,ZETA,3,1,m));
    row3.push_back(calcularTenedor(e,ZETA,4,1,m));

	matriz.push_back(row1); matriz.push_back(row2); matriz.push_back(row3);

    return determinant(matriz);
}

Matrix createLocalK(int e,mesh &m){
    
    //Preparaci�n de ingredientes -> Listo
    float u_bar,nu,rho,Ve,J,D;
    Matrix alpha,beta,gamma,delta;
    Matrix K,e_star,e_star_t,A,B,At,Bt,C,Ct;

    //Preparando alpha -> Listo
    u_bar = m.getParameter(ADJECTIVE_VELOCITY);
    J = calculateLocalJ(e,m);
    D = calculateLocalD(e,m);
    
    if(D == 0){
        cout << "\n!---CATASTROPHIC FAILURE---!\n";
        exit(EXIT_FAILURE);
    }
    
    calculateEStar(e_star);
    calculateLocalA(e,A,m);
    calculateB(B);
    
    
    productRealMatrix(u_bar*J/D,productMatrixMatrix(e_star,productMatrixMatrix(A,B,3,3,12),12,3,12),alpha);

    //Preparando beta -> Listo
    nu = m.getParameter(DYNAMIC_VISCOSITY);
    Ve = calculateLocalVolumen(e,m);

    transpose(A,At);
    transpose(B,Bt);

    productRealMatrix(nu*Ve/(D*D),productMatrixMatrix(Bt,productMatrixMatrix(At,productMatrixMatrix(A,B,3,3,12),3,3,12),12,3,12),beta);

    //Preparando gamma -> Listo
    rho = m.getParameter(DENSITY);
    calculateC(C);

    productRealMatrix(J/(rho*D),productMatrixMatrix(e_star,productMatrixMatrix(A,C,3,3,4),12,3,4),gamma);

    //Preparando delta -> Listo
    transpose(C,Ct);
    transpose(e_star,e_star_t);
    
    productRealMatrix(J/D,productMatrixMatrix(Ct,productMatrixMatrix(At,e_star_t,3,3,12),4,3,12),delta);

    //Colocando submatrices en K -> Listo!
    zeroes(K,16);
    ubicarSubMatriz(K,0,11,0,11,sumMatrix(alpha,beta,12,12));
    ubicarSubMatriz(K,0,11,12,15,gamma);
    ubicarSubMatriz(K,12,15,0,11,delta);

    return K;
}

Vector createLocalb(int e,mesh &m){
    Vector b0,b,f;
    Matrix e_star;

    float f_x = m.getParameter(EXTERNAL_FORCE_X);
    float f_y = m.getParameter(EXTERNAL_FORCE_Y);
    float f_z = m.getParameter(EXTERNAL_FORCE_Z);

    float J = calculateLocalJ(e,m);
    calculateEStar(e_star);

    zeroes(f,3);
    f.at(0) = f_x;
    f.at(1) = f_y;
    f.at(2) = f_z;

    zeroes(b0,16);

    productMatrixVector(e_star,f,b0);
    productRealVector(J,b0,b);
    
    b.push_back(0); b.push_back(0); b.push_back(0);

    return b;
}

void crearSistemasLocales(mesh &m,vector<Matrix> &localKs,vector<Vector> &localbs){
    for(int i=0;i<m.getSize(ELEMENTS);i++){
        localKs.push_back(createLocalK(i,m));
        localbs.push_back(createLocalb(i,m));
    }
}

void assemblyK(element e,Matrix localK,Matrix &K,int nnodes){
    int index1 = e.getNode1() - 1;
    int index2 = e.getNode2() - 1;
    int index3 = e.getNode3() - 1;
    int index4 = e.getNode4() - 1;
    
    int index5 = index1+nnodes;
    int index6 = index2+nnodes;
    int index7 = index3+nnodes;
    int index8 = index4+nnodes;

    int index9 = index1+2*nnodes;
    int index10 = index2+2*nnodes;
    int index11 = index3+2*nnodes;
    int index12 = index4+2*nnodes;

    int index13 = index1+3*nnodes;
    int index14 = index2+3*nnodes;
    int index15 = index3+3*nnodes;
    int index16 = index4+3*nnodes;

    K.at(index1).at(index1)  += localK.at(0).at(0);
    K.at(index1).at(index2)  += localK.at(0).at(1);
    K.at(index1).at(index3)  += localK.at(0).at(2);
    K.at(index1).at(index4)  += localK.at(0).at(3);
    K.at(index1).at(index5)  += localK.at(0).at(4);
    K.at(index1).at(index6)  += localK.at(0).at(5);
    K.at(index1).at(index7)  += localK.at(0).at(6);
    K.at(index1).at(index8)  += localK.at(0).at(7);
    K.at(index1).at(index9)  += localK.at(0).at(8);
    K.at(index1).at(index10) += localK.at(0).at(9);
    K.at(index1).at(index11) += localK.at(0).at(10);
    K.at(index1).at(index12) += localK.at(0).at(11);
    K.at(index1).at(index13) += localK.at(0).at(12);
    K.at(index1).at(index14) += localK.at(0).at(13);
    K.at(index1).at(index15) += localK.at(0).at(14);
    K.at(index1).at(index16) += localK.at(0).at(15);

    K.at(index2).at(index1)  += localK.at(1).at(0);
    K.at(index2).at(index2)  += localK.at(1).at(1);
    K.at(index2).at(index3)  += localK.at(1).at(2);
    K.at(index2).at(index4)  += localK.at(1).at(3);
    K.at(index2).at(index5)  += localK.at(1).at(4);
    K.at(index2).at(index6)  += localK.at(1).at(5);
    K.at(index2).at(index7)  += localK.at(1).at(6);
    K.at(index2).at(index8)  += localK.at(1).at(7);
    K.at(index2).at(index9)  += localK.at(1).at(8);
    K.at(index2).at(index10) += localK.at(1).at(9);
    K.at(index2).at(index11) += localK.at(1).at(10);
    K.at(index2).at(index12) += localK.at(1).at(11);
    K.at(index2).at(index13) += localK.at(1).at(12);
    K.at(index2).at(index14) += localK.at(1).at(13);
    K.at(index2).at(index15) += localK.at(1).at(14);
    K.at(index2).at(index16) += localK.at(1).at(15);
    
    K.at(index3).at(index1)  += localK.at(2).at(0);
    K.at(index3).at(index2)  += localK.at(2).at(1);
    K.at(index3).at(index3)  += localK.at(2).at(2);
    K.at(index3).at(index4)  += localK.at(2).at(3);
    K.at(index3).at(index5)  += localK.at(2).at(4);
    K.at(index3).at(index6)  += localK.at(2).at(5);
    K.at(index3).at(index7)  += localK.at(2).at(6);
    K.at(index3).at(index8)  += localK.at(2).at(7);
    K.at(index3).at(index9)  += localK.at(2).at(8);
    K.at(index3).at(index10) += localK.at(2).at(9);
    K.at(index3).at(index11) += localK.at(2).at(10);
    K.at(index3).at(index12) += localK.at(2).at(11);
    K.at(index3).at(index13) += localK.at(2).at(12);
    K.at(index3).at(index14) += localK.at(2).at(13);
    K.at(index3).at(index15) += localK.at(2).at(14);
    K.at(index3).at(index16) += localK.at(2).at(15);

    K.at(index4).at(index1)  += localK.at(3).at(0);
    K.at(index4).at(index2)  += localK.at(3).at(1);
    K.at(index4).at(index3)  += localK.at(3).at(2);
    K.at(index4).at(index4)  += localK.at(3).at(3);
    K.at(index4).at(index5)  += localK.at(3).at(4);
    K.at(index4).at(index6)  += localK.at(3).at(5);
    K.at(index4).at(index7)  += localK.at(3).at(6);
    K.at(index4).at(index8)  += localK.at(3).at(7);
    K.at(index4).at(index9)  += localK.at(3).at(8);
    K.at(index4).at(index10) += localK.at(3).at(9);
    K.at(index4).at(index11) += localK.at(3).at(10);
    K.at(index4).at(index12) += localK.at(3).at(11);
    K.at(index4).at(index13) += localK.at(3).at(12);
    K.at(index4).at(index14) += localK.at(3).at(13);
    K.at(index4).at(index15) += localK.at(3).at(14);
    K.at(index4).at(index16) += localK.at(3).at(15);

    K.at(index5).at(index1)  += localK.at(4).at(0);
    K.at(index5).at(index2)  += localK.at(4).at(1);
    K.at(index5).at(index3)  += localK.at(4).at(2);
    K.at(index5).at(index4)  += localK.at(4).at(3);
    K.at(index5).at(index5)  += localK.at(4).at(4);
    K.at(index5).at(index6)  += localK.at(4).at(5);
    K.at(index5).at(index7)  += localK.at(4).at(6);
    K.at(index5).at(index8)  += localK.at(4).at(7);
    K.at(index5).at(index9)  += localK.at(4).at(8);
    K.at(index5).at(index10) += localK.at(4).at(9);
    K.at(index5).at(index11) += localK.at(4).at(10);
    K.at(index5).at(index12) += localK.at(4).at(11);
    K.at(index5).at(index13) += localK.at(4).at(12);
    K.at(index5).at(index14) += localK.at(4).at(13);
    K.at(index5).at(index15) += localK.at(4).at(14);
    K.at(index5).at(index16) += localK.at(4).at(15);

    K.at(index6).at(index1)  += localK.at(5).at(0);
    K.at(index6).at(index2)  += localK.at(5).at(1);
    K.at(index6).at(index3)  += localK.at(5).at(2);
    K.at(index6).at(index4)  += localK.at(5).at(3);
    K.at(index6).at(index5)  += localK.at(5).at(4);
    K.at(index6).at(index6)  += localK.at(5).at(5);
    K.at(index6).at(index7)  += localK.at(5).at(6);
    K.at(index6).at(index8)  += localK.at(5).at(7);
    K.at(index6).at(index9)  += localK.at(5).at(8);
    K.at(index6).at(index10) += localK.at(5).at(9);
    K.at(index6).at(index11) += localK.at(5).at(10);
    K.at(index6).at(index12) += localK.at(5).at(11);
    K.at(index6).at(index13) += localK.at(5).at(12);
    K.at(index6).at(index14) += localK.at(5).at(13);
    K.at(index6).at(index15) += localK.at(5).at(14);
    K.at(index6).at(index16) += localK.at(5).at(15);

	K.at(index7).at(index1)  += localK.at(6).at(0);
    K.at(index7).at(index2)  += localK.at(6).at(1);
    K.at(index7).at(index3)  += localK.at(6).at(2);
    K.at(index7).at(index4)  += localK.at(6).at(3);
    K.at(index7).at(index5)  += localK.at(6).at(4);
    K.at(index7).at(index6)  += localK.at(6).at(5);
    K.at(index7).at(index7)  += localK.at(6).at(6);
    K.at(index7).at(index8)  += localK.at(6).at(7);
    K.at(index7).at(index9)  += localK.at(6).at(8);
    K.at(index7).at(index10) += localK.at(6).at(9);
    K.at(index7).at(index11) += localK.at(6).at(10);
    K.at(index7).at(index12) += localK.at(6).at(11);
    K.at(index7).at(index13) += localK.at(6).at(12);
    K.at(index7).at(index14) += localK.at(6).at(13);
    K.at(index7).at(index15) += localK.at(6).at(14);
    K.at(index7).at(index16) += localK.at(6).at(15);

	K.at(index8).at(index1)  += localK.at(7).at(0);
    K.at(index8).at(index2)  += localK.at(7).at(1);
    K.at(index8).at(index3)  += localK.at(7).at(2);
    K.at(index8).at(index4)  += localK.at(7).at(3);
    K.at(index8).at(index5)  += localK.at(7).at(4);
    K.at(index8).at(index6)  += localK.at(7).at(5);
    K.at(index8).at(index7)  += localK.at(7).at(6);
    K.at(index8).at(index8)  += localK.at(7).at(7);
    K.at(index8).at(index9)  += localK.at(7).at(8);
    K.at(index8).at(index10) += localK.at(7).at(9);
    K.at(index8).at(index11) += localK.at(7).at(10);
    K.at(index8).at(index12) += localK.at(7).at(11);
    K.at(index8).at(index13) += localK.at(7).at(12);
    K.at(index8).at(index14) += localK.at(7).at(13);
    K.at(index8).at(index15) += localK.at(7).at(14);
    K.at(index8).at(index16) += localK.at(7).at(15);

	K.at(index9).at(index1)  += localK.at(8).at(0);
    K.at(index9).at(index2)  += localK.at(8).at(1);
    K.at(index9).at(index3)  += localK.at(8).at(2);
    K.at(index9).at(index4)  += localK.at(8).at(3);
    K.at(index9).at(index5)  += localK.at(8).at(4);
    K.at(index9).at(index6)  += localK.at(8).at(5);
    K.at(index9).at(index7)  += localK.at(8).at(6);
    K.at(index9).at(index8)  += localK.at(8).at(7);
    K.at(index9).at(index9)  += localK.at(8).at(8);
    K.at(index9).at(index10) += localK.at(8).at(9);
    K.at(index9).at(index11) += localK.at(8).at(10);
    K.at(index9).at(index12) += localK.at(8).at(11);
    K.at(index9).at(index13) += localK.at(8).at(12);
    K.at(index9).at(index14) += localK.at(8).at(13);
    K.at(index9).at(index15) += localK.at(8).at(14);
    K.at(index9).at(index16) += localK.at(8).at(15);

	K.at(index10).at(index1)  += localK.at(9).at(0);
    K.at(index10).at(index2)  += localK.at(9).at(1);
    K.at(index10).at(index3)  += localK.at(9).at(2);
    K.at(index10).at(index4)  += localK.at(9).at(3);
    K.at(index10).at(index5)  += localK.at(9).at(4);
    K.at(index10).at(index6)  += localK.at(9).at(5);
    K.at(index10).at(index7)  += localK.at(9).at(6);
    K.at(index10).at(index8)  += localK.at(9).at(7);
    K.at(index10).at(index9)  += localK.at(9).at(8);
    K.at(index10).at(index10) += localK.at(9).at(9);
    K.at(index10).at(index11) += localK.at(9).at(10);
    K.at(index10).at(index12) += localK.at(9).at(11);
    K.at(index10).at(index13) += localK.at(9).at(12);
    K.at(index10).at(index14) += localK.at(9).at(13);
    K.at(index10).at(index15) += localK.at(9).at(14);
    K.at(index10).at(index16) += localK.at(9).at(15);

	K.at(index11).at(index1)  += localK.at(10).at(0);
    K.at(index11).at(index2)  += localK.at(10).at(1);
    K.at(index11).at(index3)  += localK.at(10).at(2);
    K.at(index11).at(index4)  += localK.at(10).at(3);
    K.at(index11).at(index5)  += localK.at(10).at(4);
    K.at(index11).at(index6)  += localK.at(10).at(5);
    K.at(index11).at(index7)  += localK.at(10).at(6);
    K.at(index11).at(index8)  += localK.at(10).at(7);
    K.at(index11).at(index9)  += localK.at(10).at(8);
    K.at(index11).at(index10) += localK.at(10).at(9);
    K.at(index11).at(index11) += localK.at(10).at(10);
    K.at(index11).at(index12) += localK.at(10).at(11);
    K.at(index11).at(index13) += localK.at(10).at(12);
    K.at(index11).at(index14) += localK.at(10).at(13);
    K.at(index11).at(index15) += localK.at(10).at(14);
    K.at(index11).at(index16) += localK.at(10).at(15);

	K.at(index12).at(index1)  += localK.at(11).at(0);
    K.at(index12).at(index2)  += localK.at(11).at(1);
    K.at(index12).at(index3)  += localK.at(11).at(2);
    K.at(index12).at(index4)  += localK.at(11).at(3);
    K.at(index12).at(index5)  += localK.at(11).at(4);
    K.at(index12).at(index6)  += localK.at(11).at(5);
    K.at(index12).at(index7)  += localK.at(11).at(6);
    K.at(index12).at(index8)  += localK.at(11).at(7);
    K.at(index12).at(index9)  += localK.at(11).at(8);
    K.at(index12).at(index10) += localK.at(11).at(9);
    K.at(index12).at(index11) += localK.at(11).at(10);
    K.at(index12).at(index12) += localK.at(11).at(11);
    K.at(index12).at(index13) += localK.at(11).at(12);
    K.at(index12).at(index14) += localK.at(11).at(13);
    K.at(index12).at(index15) += localK.at(11).at(14);
    K.at(index12).at(index16) += localK.at(11).at(15);

	K.at(index13).at(index1)  += localK.at(12).at(0);
    K.at(index13).at(index2)  += localK.at(12).at(1);
    K.at(index13).at(index3)  += localK.at(12).at(2);
    K.at(index13).at(index4)  += localK.at(12).at(3);
    K.at(index13).at(index5)  += localK.at(12).at(4);
    K.at(index13).at(index6)  += localK.at(12).at(5);
    K.at(index13).at(index7)  += localK.at(12).at(6);
    K.at(index13).at(index8)  += localK.at(12).at(7);
    K.at(index13).at(index9)  += localK.at(12).at(8);
    K.at(index13).at(index10) += localK.at(12).at(9);
    K.at(index13).at(index11) += localK.at(12).at(10);
    K.at(index13).at(index12) += localK.at(12).at(11);
    K.at(index13).at(index13) += localK.at(12).at(12);
    K.at(index13).at(index14) += localK.at(12).at(13);
    K.at(index13).at(index15) += localK.at(12).at(14);
    K.at(index13).at(index16) += localK.at(12).at(15);

	K.at(index14).at(index1)  += localK.at(13).at(0);
    K.at(index14).at(index2)  += localK.at(13).at(1);
    K.at(index14).at(index3)  += localK.at(13).at(2);
    K.at(index14).at(index4)  += localK.at(13).at(3);
    K.at(index14).at(index5)  += localK.at(13).at(4);
    K.at(index14).at(index6)  += localK.at(13).at(5);
    K.at(index14).at(index7)  += localK.at(13).at(6);
    K.at(index14).at(index8)  += localK.at(13).at(7);
    K.at(index14).at(index9)  += localK.at(13).at(8);
    K.at(index14).at(index10) += localK.at(13).at(9);
    K.at(index14).at(index11) += localK.at(13).at(10);
    K.at(index14).at(index12) += localK.at(13).at(11);
    K.at(index14).at(index13) += localK.at(13).at(12);
    K.at(index14).at(index14) += localK.at(13).at(13);
    K.at(index14).at(index15) += localK.at(13).at(14);
    K.at(index14).at(index16) += localK.at(13).at(15);

	K.at(index15).at(index1)  += localK.at(14).at(0);
    K.at(index15).at(index2)  += localK.at(14).at(1);
    K.at(index15).at(index3)  += localK.at(14).at(2);
    K.at(index15).at(index4)  += localK.at(14).at(3);
    K.at(index15).at(index5)  += localK.at(14).at(4);
    K.at(index15).at(index6)  += localK.at(14).at(5);
    K.at(index15).at(index7)  += localK.at(14).at(6);
    K.at(index15).at(index8)  += localK.at(14).at(7);
    K.at(index15).at(index9)  += localK.at(14).at(8);
    K.at(index15).at(index10) += localK.at(14).at(9);
    K.at(index15).at(index11) += localK.at(14).at(10);
    K.at(index15).at(index12) += localK.at(14).at(11);
    K.at(index15).at(index13) += localK.at(14).at(12);
    K.at(index15).at(index14) += localK.at(14).at(13);
    K.at(index15).at(index15) += localK.at(14).at(14);
    K.at(index15).at(index16) += localK.at(14).at(15);

	K.at(index16).at(index1)  += localK.at(15).at(0);
    K.at(index16).at(index2)  += localK.at(15).at(1);
    K.at(index16).at(index3)  += localK.at(15).at(2);
    K.at(index16).at(index4)  += localK.at(15).at(3);
    K.at(index16).at(index5)  += localK.at(15).at(4);
    K.at(index16).at(index6)  += localK.at(15).at(5);
    K.at(index16).at(index7)  += localK.at(15).at(6);
    K.at(index16).at(index8)  += localK.at(15).at(7);
    K.at(index16).at(index9)  += localK.at(15).at(8);
    K.at(index16).at(index10) += localK.at(15).at(9);
    K.at(index16).at(index11) += localK.at(15).at(10);
    K.at(index16).at(index12) += localK.at(15).at(11);
    K.at(index16).at(index13) += localK.at(15).at(12);
    K.at(index16).at(index14) += localK.at(15).at(13);
    K.at(index16).at(index15) += localK.at(15).at(14);
    K.at(index16).at(index16) += localK.at(15).at(15);

}

void assemblyb(element e,Vector localb,Vector &b,int nnodes){
    int index1 = e.getNode1() - 1;
    int index2 = e.getNode2() - 1;
    int index3 = e.getNode3() - 1;
    int index4 = e.getNode4() - 1;
    
    int index5 = index1+nnodes;
    int index6 = index2+nnodes;
    int index7 = index3+nnodes;
    int index8 = index4+nnodes;

    int index9 = index1+2*nnodes;
    int index10 = index2+2*nnodes;
    int index11 = index3+2*nnodes;
    int index12 = index4+2*nnodes;

    b.at(index1) += localb.at(0);
    b.at(index2) += localb.at(1);
    b.at(index3) += localb.at(2);
    b.at(index4) += localb.at(3);
    b.at(index5) += localb.at(4);
    b.at(index6) += localb.at(5);
    b.at(index7) += localb.at(6);
    b.at(index8) += localb.at(7);
    b.at(index9) += localb.at(8);
    b.at(index10) += localb.at(9);
    b.at(index11) += localb.at(10);
    b.at(index12) += localb.at(11);
}

void ensamblaje(mesh &m,vector<Matrix> &localKs,vector<Vector> &localbs,Matrix &K,Vector &b){
    int nnodes = m.getSize(NODES);
    for(int i=0;i<m.getSize(ELEMENTS);i++){
        element e = m.getElement(i);
        assemblyK(e,localKs.at(i),K,nnodes);
        assemblyb(e,localbs.at(i),b,nnodes);
    }
}

void applyDirichlet(mesh &m,Matrix &K,Vector &b){
    for(int i=0;i<m.getSize(DIRICHLET);i++){
        condition c = m.getCondition(i,DIRICHLET);
        int index = c.getNode1()-1;

        K.erase(K.begin()+index);
        b.erase(b.begin()+index);

        for(int row=0;row<K.size();row++){
            float cell = K.at(row).at(index);
            K.at(row).erase(K.at(row).begin()+index);
            b.at(row) += -1*c.getValue()*cell;
        }
    }
}

void calculate(Matrix &K, Vector &b, vector<Vector> &Ts,int size){
    Vector v1,v2,T,dT,Tnext;

    zeroes(v1,size);
    zeroes(v2,size);
    zeroes(T,size);
    zeroes(dT,size);

    float tf = 3;
    float delta_t = 0.5;
    float t = 0.5;
    Ts.push_back(T);

    do{
        //Tnext=T+dt*(b-K*T) //Forward Euler

        productMatrixVector(K,T,v1);
        productRealVector(-1,v1,v2);
        productRealVector(delta_t,sumVector(b,v2,size),dT);

        Tnext.clear();
        copyVector(sumVector(T,dT,size),Tnext);
        Ts.push_back(Tnext);

        T.clear();
        copyVector(Tnext,T);

        t+=delta_t;
		
    }while(t<=tf);
}