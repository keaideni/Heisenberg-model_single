#ifndef SUB_H
#define SUB_H
#include "OP.h"

class Sub
{
private:
        int _Orbital;

        OP _SubSys;
        OP _SubSysEye;
        OP _SubSysSpinP;
        OP _SubSysSpinM;
        OP _SubSysSpinZ;
        OP _SubSysSpinP1;
        OP _SubSysSpinM1;
        OP _SubSysSpinZ1;

public:
        const int& Orbital()const{return _Orbital;};
        const OP& SubSys()const{return _SubSys;};
        const OP& SubSysEye()const{return _SubSysEye;};
        const OP& SubSysSpinP()const{return _SubSysSpinP;};
        const OP& SubSysSpinM()const{return _SubSysSpinM;};
        const OP& SubSysSpinZ()const{return _SubSysSpinZ;};
        const OP& SubSysSpinP1()const{return _SubSysSpinP1;};
        const OP& SubSysSpinM1()const{return _SubSysSpinM1;};
        const OP& SubSysSpinZ1()const{return _SubSysSpinZ1;};


        Sub(){};
        ~Sub(){};


        Sub(const Sub& block);
        Sub(const int& orbital_);
        void Initial(const int& orbital);
        //=====================================================================================
        //update the whole Sub. the constant coup is for the couple constant.
        Sub(const int& orbital_, const Sub& oldL, const Sub& oldR, const double& coup);
        //======================================================================================


        void update(const int& orbital_, const Sub& oldL, const Sub& oldR, const double& coup);
        //void update(const Parameter& para, const Sub& old, const int orbital_, const int& way);//for the periodic condition.
        //truncated U to update the Sub.
        void trunc(const OP& truncU);


        void updateBlockH(const Sub& oldL, const Sub& oldR, const double& coup);//update the SubSys
        //void updateBlockH(OP& NewSys, const Parameter& para, const Sub& old, const int& way);//for the periodic condition.



        void operator=(const Sub& block);



        void show() const;
        void clear();
        void save()const;
        void read(int orbital_);


        
};


#endif // SUB_H
