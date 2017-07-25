#include "Super.h"
#include <ctime>





Super::~Super()
{

}






Super::Super(const Parameter& para_, const Sub& sys, const Sub& m, const Sub& n, const Sub& env, const int& TotQ):
RealSysBlock(sys),
RealMBlock(m),
RealNBlock(n),
RealEnvBlock(env)
{
        Wave.Initial(sys.SubSysEye(), m.SubSysEye(), n.SubSysEye(), env.SubSysEye(), TotQ);
        SaveWave = Wave;

        para = para_;

        Dim = Wave.getDim();

        
}



/*void Super::Initial(const Parameter& para_, const Sub& sys, const Sub& m, const Sub& n, const Sub& env, int TotQ)
{
        Wave = QWave(sys.SubSysEye(), m.SubSysEye(), n.SubSysEye(), env.SubSysEye(), TotQ);
        SaveWave = Wave;

        para = para_;

        Dim = Wave.getDim();

        RealSysBlock = &sys;
        RealMBlock = &m;
        RealNBlock = &n;
        RealEnvBlock = &env;
}*/



//the point f translate to g
void Super::f1tof2(const double* f, double* g)
{

        std::vector<double> f1, g1;
        for (int i = 0; i<Dim; i++)
        {
                f1.push_back(f[i]);
        }

        f1tof2(f1, g1);

        for (int i = 0; i<Dim; i++)
        {
                g[i] = g1[i];
        }
}



//the vector f translate to g.
void Super::f1tof2(const std::vector<double>& f, std::vector<double>& g)
{

        g.clear();
        Wave.f2Wave(f);
        SaveWave.setZero();
        OneStep(); //for the |SaveWave> = H |Wave>
        SaveWave.Wave2f(g);
}





//===========to calculate the QWave after the operation of the hamiltion=============
void Super::OneStep()
{

        OP p1, p2;
        //=========the term have no relationship to the boundary condition========
        Wave.OPWave(SaveWave, RealSysBlock.SubSys(), 1);//sys

        Wave.OPWave(SaveWave, RealMBlock.SubSys(), 2);//M

        Wave.OPWave(SaveWave, RealNBlock.SubSys(), 3);//N

        Wave.OPWave(SaveWave, RealEnvBlock.SubSys(), 4);//env




        //Sys-M
        p1.time(0.5, RealSysBlock.SubSysSpinM());
        p2.time(0.5, RealSysBlock.SubSysSpinP());
        BlockWave(1, p1, 2, RealMBlock.SubSysSpinP());
        BlockWave(1, p2, 2, RealMBlock.SubSysSpinM());

        BlockWave(1, RealSysBlock.SubSysSpinZ(), 2, RealMBlock.SubSysSpinZ());



                
        //==================periodic condition=====================
                
        //M-Env
        p1.time(0.5, RealMBlock.SubSysSpinM());
        p2.time(0.5, RealMBlock.SubSysSpinP());

        BlockWave(2, p1, 4, RealEnvBlock.SubSysSpinP());
        BlockWave(2, p2, 4, RealEnvBlock.SubSysSpinM());

        BlockWave(2, RealMBlock.SubSysSpinZ(), 4, RealEnvBlock.SubSysSpinZ());


                        

        //Env-N
        p1.time(0.5, RealEnvBlock.SubSysSpinM1());
        p2.time(0.5, RealEnvBlock.SubSysSpinP1());

        BlockWave(4, p1, 3, RealNBlock.SubSysSpinP());
        BlockWave(4, p2, 3, RealNBlock.SubSysSpinM());

        BlockWave(4, RealEnvBlock.SubSysSpinZ1(), 3, RealNBlock.SubSysSpinZ());

        //N-Sys
        p1.time(0.5, RealNBlock.SubSysSpinM());
        p2.time(0.5, RealNBlock.SubSysSpinP());

        BlockWave(3, p1, 1, RealSysBlock.SubSysSpinP1());
        BlockWave(3, p2, 1, RealSysBlock.SubSysSpinM1());

        BlockWave(3, RealNBlock.SubSysSpinZ(), 1, RealSysBlock.SubSysSpinZ1());
                        
                

        
        


        

}



void Super::BlockWave(const int& l, const OP& OPl, const int& r, const OP& OPr)
{
        tempWave.clear();

        tempWave.OPWave2New(Wave, OPr, r);
        tempWave.OPWave(SaveWave, OPl, l);

}



//============translate the base state to Qwave form==============
void Super::normalizedCopy(const double* f0)
{
        std::vector<double> v;
        double norm(0);

        for (int i = 0; i<Dim; i++)
        {
                norm += f0[i] * f0[i];
        }

        norm = sqrt(norm);
        for (int i = 0; i<Dim; i++)
        {
                v.push_back(f0[i] / norm);
        }

        Wave.f2Wave(v);

}




void Super::show()
{
        std::cout << "the System part: " << std::endl;
        RealSysBlock.show();
        std::cout << "the M part: " << std::endl;
        RealMBlock.show();
        std::cout << "the N part: " << std::endl;
        RealNBlock.show();
        std::cout << "the Env part: " << std::endl;
        RealEnvBlock.show();
        std::cout << "the QWave part: " << std::endl;
        Wave.show();
        std::cout << "Dim = " << Dim << std::endl;
}