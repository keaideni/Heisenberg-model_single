#include "Sub.h"

Sub::Sub(const Sub& block):
_Orbital(block._Orbital),
_SubSys(block._SubSys),
_SubSysEye(block._SubSysEye),
_SubSysSpinP(block._SubSysSpinP),
_SubSysSpinM(block._SubSysSpinM),
_SubSysSpinZ(block._SubSysSpinZ),
_SubSysSpinP1(block._SubSysSpinP1),
_SubSysSpinM1(block._SubSysSpinM1),
_SubSysSpinZ1(block._SubSysSpinZ1)
{}

Sub::Sub(const int& orbital):
_Orbital(orbital),
_SubSysEye(SpinEye),
_SubSysSpinP(SpinP),
_SubSysSpinM(SpinM),
_SubSysSpinZ(SpinZ),
_SubSysSpinP1(SpinP),
_SubSysSpinM1(SpinM),
_SubSysSpinZ1(SpinZ)
{

}


void Sub::Initial(const int& orbital)
{
        _Orbital=orbital;
        _SubSysEye.Initial(SpinEye);
        _SubSysSpinP.Initial(SpinP);
        _SubSysSpinM.Initial(SpinM);
        _SubSysSpinZ.Initial(SpinZ);
        _SubSysSpinP1.Initial(SpinP);
        _SubSysSpinM1.Initial(SpinM);
        _SubSysSpinZ1.Initial(SpinZ);

}

Sub::Sub(const int& orbital, const Sub& oldL, const Sub& oldR, const double& coup):
_Orbital(orbital)
{
        updateBlockH(oldL, oldR, coup);
        _SubSysSpinP.kron(oldL._SubSysEye, oldR._SubSysSpinP);
        _SubSysSpinM.kron(oldL._SubSysEye, oldR._SubSysSpinM);
        _SubSysSpinZ.kron(oldL._SubSysEye, oldR._SubSysSpinZ);
        _SubSysEye.kron(oldL._SubSysEye, oldR._SubSysEye);
        _SubSysSpinP1.kron(oldL._SubSysSpinP1, oldR._SubSysEye);
        _SubSysSpinM1.kron(oldL._SubSysSpinM1, oldR._SubSysEye);
        _SubSysSpinZ1.kron(oldL._SubSysSpinZ1, oldR._SubSysEye);

}


void Sub::updateBlockH(const Sub& oldL, const Sub& oldR, const double& coup)
{

        clear();
        OP tempOP1;
        //OP tempOP2;



        _SubSys.kron(oldL._SubSysEye, oldR._SubSys);
        //tempOP2.time(tempOP1, para.Wc);



        //H kron I.

        //tempOP2.clear();
        tempOP1.kron(oldL._SubSys, oldR._SubSysEye);
        _SubSys.add(tempOP1); 



        tempOP1.kron(oldL._SubSysSpinP, oldR._SubSysSpinM1);
        //oldL.SubSysCdag.show();
        //oldR.SubSysC1.show();
        //tempOP1.show();
        tempOP1.time(coup);
        _SubSys.add(tempOP1);//_SubSys.show();


        tempOP1.clear();
        tempOP1.kron(oldL._SubSysSpinM, oldR._SubSysSpinP1);
        tempOP1.time(coup);
        _SubSys.add(tempOP1);

        tempOP1.clear();
        tempOP1.kron(oldL._SubSysSpinZ, oldR._SubSysSpinZ1);
        _SubSys.add(tempOP1);




}


void Sub::update(const int& orbital, const Sub& oldL, const Sub& oldR, const double& coup)
{
        clear();
        _Orbital=orbital;
        updateBlockH(oldL, oldR, coup);
        _SubSysSpinP.kron(oldL._SubSysEye, oldR._SubSysSpinP);
        _SubSysSpinM.kron(oldL._SubSysEye, oldR._SubSysSpinM);
        _SubSysSpinZ.kron(oldL._SubSysEye, oldR._SubSysSpinZ);
        _SubSysEye.kron(oldL._SubSysEye, oldR._SubSysEye);
        _SubSysSpinP1.kron(oldL._SubSysSpinP1, oldR._SubSysEye);
        _SubSysSpinM1.kron(oldL._SubSysSpinM1, oldR._SubSysEye);
        _SubSysSpinZ1.kron(oldL._SubSysSpinZ1, oldR._SubSysEye);
}


void Sub::trunc(const OP& truncU)
{
        _SubSys.trunc(truncU);
        _SubSysEye.trunc(truncU);
        _SubSysSpinP.trunc(truncU);
        _SubSysSpinM.trunc(truncU);
        _SubSysSpinZ.trunc(truncU);
        _SubSysSpinP1.trunc(truncU);
        _SubSysSpinM1.trunc(truncU);
        _SubSysSpinZ1.trunc(truncU);

}


void Sub::operator=(const Sub& block)
{
        clear();

        _Orbital = block._Orbital;
        _SubSys = block._SubSys;
        _SubSysEye = block._SubSysEye;
        _SubSysSpinP = block._SubSysSpinP;
        _SubSysSpinM = block._SubSysSpinM;
        _SubSysSpinZ = block._SubSysSpinZ;
        _SubSysSpinP1 = block._SubSysSpinP1;
        _SubSysSpinM1 = block._SubSysSpinM1;
        _SubSysSpinZ1 = block._SubSysSpinZ1;

        


}




void Sub::show() const
{
        cout<<"===============The "<<_Orbital<<"th SubLattice:==============="<<endl;
        cout<<"==========The SubSys:=========="<<endl;
        _SubSys.show();
        cout<<"==========The SubSysSpinP:=========="<<endl;
        _SubSysSpinP.show();
        cout<<"==========The SubSysSpinM:=========="<<endl;
        _SubSysSpinM.show();
        cout<<"==========The SubSysSpinZ:=========="<<endl;
        _SubSysSpinZ.show();
        cout<<"==========The SubSysSpinP1:=========="<<endl;
        _SubSysSpinP1.show();
        cout<<"==========The SubSysSpinM1:=========="<<endl;
        _SubSysSpinM1.show();
        cout<<"==========The SubSysSpinZ1:=========="<<endl;
        _SubSysSpinZ1.show();
        cout<<"===============The end==============="<<endl;
}

void Sub::clear()
{
        _SubSys.clear();
        _SubSysEye.clear();
        _SubSysSpinP.clear();
        _SubSysSpinM.clear();
        _SubSysSpinZ.clear();
        _SubSysSpinP1.clear();
        _SubSysSpinM1.clear();
        _SubSysSpinZ1.clear();
}


void Sub::save()const
{
        std::string str = itos(_Orbital);
        str = "./data/" + str;
        std::ofstream outfile(str);
        outfile << _Orbital << std::endl;
        _SubSys.save(outfile);
        _SubSysEye.save(outfile);
        _SubSysSpinP.save(outfile);
        _SubSysSpinM.save(outfile);
        _SubSysSpinZ.save(outfile);
        _SubSysSpinP1.save(outfile);
        _SubSysSpinM1.save(outfile);
        _SubSysSpinZ1.save(outfile);
        outfile.close();

}

void Sub::read(int orbital_)
{
        clear();
        std::string str = itos(orbital_);
        str = "./data/" + str;
        std::ifstream infile(str);
        if (!infile.is_open())
        {
                std::cout << "the file " << str << " could not open. " << std::endl;
                exit(1);
        }
        infile >> _Orbital;
        _SubSys.read(infile);
        _SubSysEye.read(infile);
        _SubSysSpinP.read(infile);
        _SubSysSpinM.read(infile);
        _SubSysSpinZ.read(infile);
        _SubSysSpinP1.read(infile);
        _SubSysSpinM1.read(infile);
        _SubSysSpinZ1.read(infile);

}
