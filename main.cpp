#include "test.h"
#include "DMRG.h"

int OP::Max;


int main()
{
        Parameter para;//(40, 40, 200, 6);
        para.read();

        OP::Max=para.ParticleNo()+1;

        DMRG Heisenberg(para);
}