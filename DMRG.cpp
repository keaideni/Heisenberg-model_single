#include "DMRG.h"




DMRG::DMRG()
{

}
DMRG::~DMRG()
{

}

DMRG::DMRG(Parameter& para)
{
        clock_t Abegin;
        double allT;
        int OS(1);
        int OE(para.LatticeSize());
        int dir(1);


        //===========this label for the way of growth the blocks======
        OrbitalM = para.LatticeSize()/4+1;
        OrbitalN = para.LatticeSize()/4+1;
        Gdir = 1;
        //====================================================


        SaveAll.open("./result/SaveAll");
        if (!SaveAll.is_open())
        {
                //std::cout << "the file doesn't exit!" << std::endl;
        }

        saveT = 0;

        SaveAll << "==============build up================" << std::endl;
        Abegin = clock();





        BuildUpP(para, OS, OE, dir);

        allT = difftime(clock(), Abegin) / CLOCKS_PER_SEC;
        SaveAll << "===========build up finished============" << std::endl;

        SaveAll << "The build up process takes " << allT << " seconds!" << std::endl;
        SaveAll << "The eigenstate calculation takes " << saveT << " seconds!" << std::endl;

        SaveAll << "===============Sweep================" << std::endl;
        saveT = 0;
        Abegin = clock();


        Fdata.open("./result/data");
        if (!Fdata.is_open())
        {
                std::cout << "the file doesn't exit!" << std::endl;
        }

        calnonestepSM=calnonestepSN=calnonestepEM=calnonestepEN=0;
        calntwostepSM=calntwostepSN=calntwostepEM=calntwostepEN=1;
        caln = 0;
        SweepP(para, OS, OE, dir);

//======================================================================================================
        allT = difftime(clock(), Abegin) / CLOCKS_PER_SEC;
        SaveAll << "===========Sweep finished==============" << std::endl;
        SaveAll << "The build up process takes " << allT << " seconds!" << std::endl;
        SaveAll << "The eigenstate calculation takes " << saveT << " seconds!" << std::endl;
        SaveAll.close();
//===========================save the final wave==========================================
        fwave.save();
        //fwave.show();

        /*QWave tempwave;
        std::ifstream infile("./Corr/QWave");
        if(!infile.is_open())
        {
                std::cout<<"the ./Corr/QWave file doesn't open!"<<std::endl;
        }
        //std::cout<<"haha"<<std::endl;

        tempwave.read(infile);*/
        //tempwave.show();
        


        

        Fdata.close();
}






void DMRG::BuildUpP(Parameter& para, int& OS, int& OE, int& dir)
{

        befortruncateP(para, OS, OE);
        //int judge = (OE-OS ==1)? 1:0;
        int i(1);

        while (true)
        {
                SaveAll << i << std::endl;

                Sys.read(OS);//Sys.show();
                Env.read(OE);//Env.show();

                m.Initial(OrbitalM);
                n.Initial(OrbitalN);

                getEnergyP(para, dir);






                if ((OE - OS) == 3)break;


                truncUpdateP(para, OS, OE, dir);

                if (Gdir == 1)
                {
                        OrbitalM += 1;
                }
                else
                {
                        OrbitalN += 1;
                }
                Gdir *= -1;


                i++;
        }

}




void DMRG::befortruncateP(const Parameter& para, const int& OS, const int& OE) const
{
        Sub Sys1(OS);
        Sys1.save();
        //std::cout<<Sys.SubSys.Max<<std::endl;
        Sub Env1(OE);
        Env1.save();
}



void DMRG::getEnergyP(Parameter& para, int dir)
{



        int qtot = Sys.Orbital()*2;
        //std::cout<<qtot<<std::endl;
        Super Sup(para, Sys, m, n, Env, qtot);
        //std::cout<<"hehe"<<std::endl;
        begin = clock();
        SuperEnergy Supp(para, Sup);
        saveT += difftime(clock(), begin) / CLOCKS_PER_SEC;
        //std::cout<<"haha"<<std::endl;
        double trace;
        double truncerr;


        OP temp;
        OP dentemp;


        Supp.wave.Wave2OP(temp, Sys.SubSysEye(), m.SubSysEye(), n.SubSysEye(), Env.SubSysEye(), Gdir);

        dentemp.getDenS(temp);
        truncU.DengetTruncU(para, dentemp, trace, truncerr);//temp.show();

        temp.clear();
        Supp.wave.Wave2OP(temp, Sys.SubSysEye(), m.SubSysEye(), n.SubSysEye(), Env.SubSysEye(), -1 * Gdir);




        //Sup.Wave.Wave2OP(temp, Sys.SubSysEye, m.SubSysEye, n.SubSysEye, Env.SubSysEye, Gdir);
        dentemp.getDenE(temp);
        truncUR.DengetTruncU(para, dentemp, trace, truncerr);

        //truncU.show();
        /*SaveAll << "Q = " << qtot << "    WaveD = " << std::setw(4) << Sup.Dim
                << "      OS =" << std::setw(2) << Sys.Orbital() << ",  OE =" << std::setw(2) << Env.Orbital()
                << ",    E = " << std::setw(10) << std::setprecision(15) << para.Energy << ",    trace = " << std::setprecision(15) << trace
                << ",    truncerr = " << std::setprecision(15) << truncerr << std::endl;*/

        SaveAll << "Q="<<std::setw(4) << qtot << ",  WaveD=" <<std::setw(8)<< Sup.Dim
        << ",  OS="  <<std::setw(2)<<Sys.Orbital() << ",  OE=" <<std::setw(2)<< Env.Orbital()
        << ",  E=" <<std::setw(18)<< std::setprecision(15)<<para.Energy <<",  trace ="<<setw(18)
        <<std::setprecision(15)<<trace
        <<",  truncerr="<<setw(20)<< std::setprecision(15)<<truncerr<<std::endl;




        std::cout << "Q="<<std::setw(4) << qtot << ",  WaveD=" <<std::setw(8)<< Sup.Dim
        << ",  OS="  <<std::setw(2)<<Sys.Orbital() << ",  OE=" <<std::setw(2)<< Env.Orbital()
        << ",  E=" <<std::setw(18)<< std::setprecision(15)<<para.Energy <<",  trace ="<<setw(18)
        <<std::setprecision(15)<<trace
        <<",  truncerr="<<setw(20)<< std::setprecision(15)<<truncerr<<std::endl;
        FEnergy = para.Energy;
        FTrace = trace;
        FTruncerr = truncerr;



}




void DMRG::truncUpdateP(const Parameter& para, int& OS, int& OE, int dir)
{

        OS += dir;
        OE -= dir;

        if (Gdir == 1)
        {
               
                newS.update(OS, Sys, m, 0.5);
                newE.update(OE, Env, m, 0.5);

               

        }
        else
        {
                
                newS.update(OS, n, Sys, 0.5);
                newE.update(OE, n, Env, 0.5);
                

        }





        //if(OS != 1)
        //{
        newS.trunc(truncU);
        newE.trunc(truncUR);
        //}

        newS.save();





        newE.save();

        truncU.truncsave(OS);
        truncUR.truncsave(OE);

}



//=============sweep================
void DMRG::SweepP(Parameter& para, int& OS, int& OE, int& dir)
{

        int flag(0);
        //para.read();
        FEnergy = 10000000000;
        errEnergy=1;
        while (errEnergy>0.0001)
        {

                SaveAll << "the " << (flag + 1) << "th Sweep" << std::endl;
                std::cout<<"the "<<(flag+1)<<"th Sweep"<<std::endl;
                //dir*=(-1);//local here for the first left direction sweep


                //FEnergy = 1000000000;
                while (true)
                {
                        m.Initial(OrbitalM);
                        n.Initial(OrbitalN);
                        //===============consider the ways in the xishoudian d fangshi ==================================================================================
                        //==========saomiao guocheng zhong d zhangdian, zai zhangdao zhongdian d shihou you yige zhangdian fangxiang d fanzhuang, yuanben yinggai ============
                        //==========youbian zhangdian d shihou huancheng l zuobian zhangdian.=========================================================================
                        //=======zhuyao kan Wave2OP d fangshi====================


                        Sys.read(OS);//Sys.show();
                        Env.read(OE);//Env.show();

//====================================two step of the initial wave====================================================================
                        if(dir == 1)
                        {
                                if(Gdir == -1)
                                {
                                        if(calnonestepSM == calntwostepSM)
                                        {
                                                //ffwave.show();
                                                initwave.twostepSM(onewave, Sys.SubSysEye(), m.SubSysEye(), Env.SubSysEye(), n.SubSysEye());
                                                //initwave1.twostepSM(onewave1, Sys.SubSysEye(), m.SubSysEye(), Env.SubSysEye(), n.SubSysEye());
                                                //initwave2.twostepSM(onewave2, Sys.SubSysEye(), m.SubSysEye(), Env.SubSysEye(), n.SubSysEye());

                                                //fwave11.show();
                                                //exit(true);
                                                ++calntwostepSM;
                                        }
                                }else
                                {
                                        if(calnonestepSN == calntwostepSN)
                                        {
                                                initwave.twostepSN(onewave, Sys.SubSysEye(), m.SubSysEye(), Env.SubSysEye(), n.SubSysEye());
                                                //initwave1.twostepSN(onewave1, Sys.SubSysEye(), m.SubSysEye(), Env.SubSysEye(), n.SubSysEye());
                                                //initwave2.twostepSN(onewave2, Sys.SubSysEye(), m.SubSysEye(), Env.SubSysEye(), n.SubSysEye());

                                                //fwave11.show();exit(true);
                                                ++calntwostepSN;
                                        }
                                }

                        }else
                        {
                                if(Gdir == 1)
                                {
                                        if(calnonestepEM == calntwostepEM)
                                        {
                                                initwave.twostepEM(onewave, Sys.SubSysEye(), m.SubSysEye(), Env.SubSysEye(), n.SubSysEye());
                                                //initwave1.twostepEM(onewave1, Sys.SubSysEye(), m.SubSysEye(), Env.SubSysEye(), n.SubSysEye());
                                                //initwave2.twostepEM(onewave2, Sys.SubSysEye(), m.SubSysEye(), Env.SubSysEye(), n.SubSysEye());

                                                //fwave11.show();exit(true);
                                                ++calntwostepEM;
                                        }
                                }else
                                {
                                        if(calnonestepEN == calntwostepEN)
                                        {
                                                //ffwave.show();
                                                initwave.twostepEN(onewave, Sys.SubSysEye(), m.SubSysEye(), Env.SubSysEye(), n.SubSysEye());
                                                //initwave1.twostepEN(onewave1, Sys.SubSysEye(), m.SubSysEye(), Env.SubSysEye(), n.SubSysEye());
                                                //initwave2.twostepEN(onewave2, Sys.SubSysEye(), m.SubSysEye(), Env.SubSysEye(), n.SubSysEye());

                                                //fwave11.show();exit(true);
                                                ++calntwostepEN;
                                        }

                                }
                        }
//==================================================this position is very important================================================


                        //==================this aprt is for the first right sweep, if first left sweep, it should absent==========
                        if (OS == (para.LatticeSize()-2)/2)
                        {
                                Gdir *= -1;

                        }
                        //============================present with line dir *= -1 at the end of while(true)=================================

                        

                        


                        OP truncU;

                        //exit(1);

                        getEnergySweepP(para, dir);


                        


                        //this one is for the break point at the middle of the line.

                        if((errEnergy<0.0001) &&(OS == (para.LatticeSize() - 2) / 2))
                        {
                                




                                flag = para.SweepNo();
                                break;
                        }

                        if (dir == 1)
                        {

                                if (para.LatticeSize()==OE)
                                {

                                        break;
                                }
                        }
                        else
                        {
                                if (OS == 1)
                                {

                                        break;
                                }
                        }

                        
                        truncUpdateSweepP(para, OS, OE, dir);




                        if (dir == 1)
                        {
                                if (Gdir == 1)
                                {
                                        OrbitalM += 1;
                                }
                                else
                                {
                                        OrbitalN += 1;
                                }
                                Gdir *= -1;
                        }
                        else
                        {

                                if (Gdir == 1)
                                {
                                        OrbitalN -= 1;
                                }
                                else
                                {
                                        OrbitalM -= 1;
                                }
                                Gdir *= -1;
                        }


                        //==================this aprt is for the first left sweep, if first left sweep, it should absent==========
                        /*if(OS == (para.LatticeSize -2)/2)
                        Gdir *= -1;*/
                        //============================present with line dir *= -1 before while(true)========================================

                        ++caln;

                }
                dir *= (-1);    //local the for the first right sweep
                flag++;
        }

        Fdata << "Q = " << para.ParticleNo() << "    LatticeSize = " << std::setw(4) << para.LatticeSize() 
                << ",    E = " << std::setprecision(15) << Energy
                << ",    trace = " << std::setprecision(15) 
                << FTrace << ",    truncerr = " << std::setprecision(15) << FTruncerr 
                << "              para.D = "<<std::setprecision(15)<<para.D()
                <<"          Entanglment = "<<std::setprecision(15)<<FEntanglement<<std::endl;
        cout << "Q = " << para.ParticleNo() << "    LatticeSize = " << std::setw(4) << para.LatticeSize() 
                << ",    E = " << std::setprecision(15) << Energy
                << ",    trace = " << std::setprecision(15) 
                << FTrace << ",    truncerr = " << std::setprecision(15) << FTruncerr 
                << "              para.D = "<<std::setprecision(15)<<para.D()
                <<"          Entanglment = "<<std::setprecision(15)<<FEntanglement<<std::endl;
        
        
}



void DMRG::getEnergySweepP(Parameter& para, int dir)
{



        int qtot = para.ParticleNo();
        double trace;
        double truncerr;
        double Entanglment;




        Super Sup(para, Sys, m, n, Env, qtot);
        Super SupAdd1(para, Sys, m, n, Env, qtot+1);
        Super SupAdd2(para, Sys, m, n, Env, qtot+2);
        //Sup.Wave.show();

        //if(caln != 0)
        //{
                

        //}

        

        begin = clock();
        SuperEnergy Supp, SuppAdd1, SuppAdd2;
        if(caln == 0)//(Sys.Orbital == (para.ParticleNo-1)))
        {
                Supp.init(para, Sup);Energy = para.Energy;
                //SuppAdd1.init(para, SupAdd1); Add1Energy=para.Energy;
                //SuppAdd2.init(para, SupAdd2); Add2Energy=para.Energy;
        }else
        {

        
                //std::cout<<"Sys.Orbital"<<Sys.Orbital<<std::endl;
                
                Supp.init(para, Sup, initwave);Energy = para.Energy;//exit(true); 
                //SuppAdd1.init(para, SupAdd1, initwave1); Add1Energy=para.Energy;
                //SuppAdd2.init(para, SupAdd2, initwave2); Add2Energy=para.Energy;
                 

        }
        saveT += difftime(clock(), begin) / CLOCKS_PER_SEC;
        startwave = Supp.wave;//temp now
        //startwave1 = SuppAdd1.wave;
        //startwave2 = SuppAdd2.wave;


//====================================================================================================
        OP temp;
        //OP Dentemp;
        Supp.wave.Wave2OP(temp, Sys.SubSysEye(), m.SubSysEye(), n.SubSysEye(), Env.SubSysEye(), Gdir);
        if (dir == 1)
        {
                DenOPWave.getDenS(temp);
        }
        else
        {
                DenOPWave.getDenE(temp);
        }
        /*DenOPWave.time(1.0/3);


        temp.clear();OP tempDen;
        SuppAdd1.wave.Wave2OP(temp, Sys.SubSysEye(), m.SubSysEye(), n.SubSysEye(), Env.SubSysEye(), Gdir);
        if (dir == 1)
        {
                tempDen.getDenS(temp);
        }
        else
        {
                tempDen.getDenE(temp);
        }

        tempDen.time(1.0/3);
        DenOPWave.add(tempDen);


        temp.clear();tempDen.clear();
        //OP Dentemp;
        SuppAdd2.wave.Wave2OP(temp, Sys.SubSysEye(), m.SubSysEye(), n.SubSysEye(), Env.SubSysEye(), Gdir);
        if (dir == 1)
        {
                tempDen.getDenS(temp);
        }
        else
        {
                tempDen.getDenE(temp);
        }
        tempDen.time(1.0/3);
        DenOPWave.add(tempDen);*/
//=====================================================================================================       

        
        



        truncU.DengetTruncU(para, DenOPWave, trace, truncerr, Entanglment);//temp.show();






        //truncU.show();
        /*SaveAll << "Q = " << qtot << ",    E = " << std::setprecision(15) << Energy  << std::endl
                << "      OS =" << std::setw(2) << Sys.Orbital() << ",  OE =" << std::setw(2) << Env.Orbital()
                << ",    trace = " << std::setprecision(15) << trace
                << ",    truncerr = " << std::setprecision(15) << truncerr << std::endl << std::endl;*/

        SaveAll << "Q=" << qtot << ",  E = " <<setw(18)<< std::setprecision(15)<<Energy <<std::endl
        << "  OS="  <<std::setw(2)<<Sys.Orbital() << ",  OE=" <<std::setw(2)<< Env.Orbital()
        << "  OrbitalM ="  <<std::setw(2)<<OrbitalM << ",  OrbitalN=" <<std::setw(2)<< OrbitalN
        <<",  trace="<<setw(18)<< std::setprecision(15)<<trace
        <<",  truncerr="<<setw(20)<< std::setprecision(15)<<truncerr << std::endl<<std::endl;


        std::cout << "Q=" << qtot << ",  E = " <<setw(18)<< std::setprecision(15)<<Energy <<std::endl
        << "  OS="  <<std::setw(2)<<Sys.Orbital() << ",  OE=" <<std::setw(2)<< Env.Orbital()
        << "  OrbitalM ="  <<std::setw(2)<<OrbitalM << ",  OrbitalN=" <<std::setw(2)<< OrbitalN
        <<",  trace="<<setw(18)<< std::setprecision(15)<<trace
        <<",  truncerr="<<setw(20)<< std::setprecision(15)<<truncerr << std::endl<<std::endl;

        if (Sys.Orbital() == (para.LatticeSize() - 2) / 2)
        {
                errEnergy=abs(Energy-FEnergy);
                FEnergy = Energy;
                //FAdd1Energy=Add1Energy;
                //FAdd2Energy=Add2Energy;
                
                FTrace = trace;
                FTruncerr = truncerr;
                FEntanglement = Entanglment;
                fwave = Supp.wave;
                
        }



}



void DMRG::truncUpdateSweepP(const Parameter& para, int& OS, int& OE, int dir)
{

        //CorrUpdate(dir, para);
        /*std::cout<<"the corr Dim: "<<std::endl;
        for(auto it = corr.CorrO.QDim.begin(); it != corr.CorrO.QDim.end(); ++it)
        {
                std::cout << it->first << "=>" << it->second<<std::endl;
        }*/

        


        //std::cout<<dir<<std::endl;
        if (dir == 1)
        {

                if (Gdir == 1)
                {
                        
                        newS.update(OS + 1, Sys, m, 0.5);
                        

//============================================================================================================================================================
                        //this part is for the wavetransform.
                        //QWave ffwave;
                        OP truncE;
                        truncE.truncread(OE);//truncE.show();
                        onewave.onestepSM(startwave, Sys.SubSysEye(), m.SubSysEye(), Env.SubSysEye(), n.SubSysEye(), truncU, truncE);
                        //onewave1.onestepSM(startwave1, Sys.SubSysEye(), m.SubSysEye(), Env.SubSysEye(), n.SubSysEye(), truncU, truncE);
                        //onewave2.onestepSM(startwave2, Sys.SubSysEye(), m.SubSysEye(), Env.SubSysEye(), n.SubSysEye(), truncU, truncE);

                        ++calnonestepSM;
                        //fwave1.show();ffwave.show();exit(true);
//============================================================================================================================================================
                }
                else
                {
                        
                        newS.update(OS + 1, n, Sys, 0.5);
                        
//============================================================================================================================================================
                        //this part is for the wavetransform.
                        //QWave ffwave;
                        OP truncE;
                        truncE.truncread(OE);//truncE.show();
                        onewave.onestepSN(startwave, Sys.SubSysEye(), m.SubSysEye(), Env.SubSysEye(), n.SubSysEye(), truncU, truncE);
                        //onewave1.onestepSN(startwave1, Sys.SubSysEye(), m.SubSysEye(), Env.SubSysEye(), n.SubSysEye(), truncU, truncE);
                        //onewave2.onestepSN(startwave2, Sys.SubSysEye(), m.SubSysEye(), Env.SubSysEye(), n.SubSysEye(), truncU, truncE);

                        ++calnonestepSN;
                        //ffwave.show();//exit(true);
//=============================================================================================================================================================
                }
                //newS.show();
                //truncU.show();

                newS.trunc(truncU);
                //newS.show();

                truncU.truncsave(OS+1);

                newS.save();
        }
        else
        {
                


                if (Gdir == 1)
                {
                        
                        newS.update(OE - 1, n, Env, 0.5);
                       
//============================================================================================================================================================
                        //this part is for the wavetransform.
                        //QWave ffwave;
                        OP truncS;
                        truncS.truncread(OS);//truncE.show();
                        onewave.onestepEN(startwave, Sys.SubSysEye(), m.SubSysEye(), Env.SubSysEye(), n.SubSysEye(), truncU, truncS);
                        //onewave1.onestepEN(startwave1, Sys.SubSysEye(), m.SubSysEye(), Env.SubSysEye(), n.SubSysEye(), truncU, truncS);
                        //onewave2.onestepEN(startwave2, Sys.SubSysEye(), m.SubSysEye(), Env.SubSysEye(), n.SubSysEye(), truncU, truncS);

                        ++calnonestepEN;
                        //ffwave.show();exit(true);
//============================================================================================================================================================
                }
                else
                {


                        
                        newS.update(OE - 1, Env, m, 0.5);
                       
//============================================================================================================================================================
                        //this part is for the wavetransform.
                        //QWave ffwave;
                        OP truncS;
                        truncS.truncread(OS);//truncE.show();
                        onewave.onestepEM(startwave, Sys.SubSysEye(), m.SubSysEye(), Env.SubSysEye(), n.SubSysEye(), truncU, truncS);
                        //onewave1.onestepEM(startwave1, Sys.SubSysEye(), m.SubSysEye(), Env.SubSysEye(), n.SubSysEye(), truncU, truncS);
                        //onewave2.onestepEM(startwave2, Sys.SubSysEye(), m.SubSysEye(), Env.SubSysEye(), n.SubSysEye(), truncU, truncS);

                        ++calnonestepEM;
                        //ffwave.show();exit(true);
//============================================================================================================================================================
                        
                }




                newS.trunc(truncU);
                //truncUR.show();
                truncU.truncsave(OE-1);

                newS.save();
        }

        /*std::cout<<"the newS Dim: "<<std::endl;
        for(auto it = newS.SubSysEye.QDim.begin(); it != newS.SubSysEye.QDim.end(); ++it)
        {
                std::cout << it->first << "=>" << it->second<<std::endl;
        }*/
        OE += dir;
        OS += dir;

        //fwave.show();
        
}