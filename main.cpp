#include <iostream>
#include "math_tools.h"
#include "classes.h"
#include "tools.h"
#include "sel.h"

int main(int argc, char *argv[])
{
    char filename[150];
    strcpy(filename,argv[1]);

    vector<Matrix> localKs;
    vector<Vector> localbs;
    vector<Vector> Vs;
    Matrix K;
    Vector b;

    cout << "IMPLEMENTACI"<<char(224)<<"N DEL M"<<char(144)<<"TODO DE LOS ELEMENTOS FINITOS\n"
         << "\t- ECUACIONES DE NAVIER-STOKES\n" << "\t- 3 DIMENSIONES\n"
         << "\t- FUNCIONES DE FORMA LINEALES\n" << "\t- PESOS DE GALERKIN\n"
         << "*********************************************************************************\n\n";

    mesh m;
    leerMallayCondiciones(m,filename);
    cout << "Datos obtenidos correctamente\n********************\n";

    crearSistemasLocales(m,localKs,localbs);
    cout<<"creo sistemas locales"<<endl;
    //showKs(localKs); showbs(localbs);
    cout << "******************************\n";

    zeroes(K,4*m.getSize(NODES));
    zeroes(b,4*m.getSize(NODES));
    
    ensamblaje(m,localKs,localbs,K,b);

    //showMatrix(K); showVector(b);
    cout << "******Ensamblaje listo***********\n";


    applyDirichlet(m,K,b);
    
    //showMatrix(K); showVector(b);
    cout << "**********Dirichlet Aplicado*********\n";

    calculate(K,b,Vs,4*m.getSize(NODES)-m.getSize(DIRICHLET));
    cout<<"calculos listos"<<endl;
    //showbs(Vs);

    writeResults(m,Vs,filename);
	cout<<"archivo creado"<<endl;
    return 0;
}
