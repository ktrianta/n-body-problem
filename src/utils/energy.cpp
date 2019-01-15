#include <energy.hpp>

sim::data_type energy(const size_t N, sim::data_type (*r)[7]) {

sim::data_type kinEn=0;
sim::data_type potEn=0;
sim::data_type en;
for (int i = 0; i < N; i++){
    kinEn += r[i][0] * (r[i][4]*r[i][4] + r[i][5]*r[i][5] + r[i][6]*r[i][6])/2.;
    for (int j = 0; j < i; j++){
        sim::data_type denominator = sqrt((r[j][1]-r[i][1])*(r[j][1]-r[i][1]) + (r[j][2]-r[i][2])*(r[j][2]-r[i][2]) +
                                  (r[j][3]-r[i][3])*(r[j][3]-r[i][3]));
        potEn -= sim::g*r[i][0]*r[j][0]/denominator;
        }
    }
    en = kinEn + potEn;

return en;
}
sim::data_type energy(const size_t N, sim::data_type (*r)[8]) {

sim::data_type kinEn=0;
sim::data_type potEn=0;
sim::data_type en;
for (int i = 0; i < N; i++){
    kinEn += r[i][0] * (r[i][4]*r[i][4] + r[i][5]*r[i][5] + r[i][6]*r[i][6])/2.;
    for (int j = 0; j < i; j++){
        sim::data_type denominator = sqrt((r[j][1]-r[i][1])*(r[j][1]-r[i][1]) + (r[j][2]-r[i][2])*(r[j][2]-r[i][2]) +
                                  (r[j][3]-r[i][3])*(r[j][3]-r[i][3]));
        potEn -= sim::g*r[i][0]*r[j][0]/denominator;
        }
    }
    en = kinEn + potEn;

return en;
}

sim::data_type energy(const size_t N, sim::data_type (*r)[3], sim::data_type (*u)[3], sim::data_type *m) {

sim::data_type kinEn=0;
sim::data_type potEn=0;
sim::data_type en;
for (int i = 0; i < N; i++){
    kinEn += m[i] * (u[i][0]*u[i][0] + u[i][1]*u[i][1] + u[i][2]*u[i][2])/2.;
    for (int j = 0; j < i; j++){
        sim::data_type denominator = sqrt((r[j][0]-r[i][0])*(r[j][0]-r[i][0]) + (r[j][1]-r[i][1])*(r[j][1]-r[i][1]) +
                                  (r[j][2]-r[i][2])*(r[j][2]-r[i][2]));
        potEn -= sim::g*m[i]*m[j]/denominator;
        }
    }
    en = kinEn + potEn;

return en;
}

sim::data_type energy(const size_t N, sim::data_type *rx, sim::data_type *ry, sim::data_type *rz, sim::data_type (*u)[3], sim::data_type *m) {

sim::data_type kinEn=0;
sim::data_type potEn=0;
sim::data_type en;
for (int i = 0; i < N; i++){
    kinEn += m[i] * (u[i][0]*u[i][0] + u[i][1]*u[i][1] + u[i][2]*u[i][2])/2.;
    for (int j = 0; j < i; j++){
        sim::data_type denominator = sqrt((rx[j]-rx[i])*(rx[j]-rx[i]) + (ry[j]-ry[i])*(ry[j]-ry[i]) +
                                  (rz[j]-rz[i])*(rz[j]-rz[i]));
        potEn -= sim::g*m[i]*m[j]/denominator;
        }
    }
    en = kinEn + potEn;

return en;
}
void printEnergy(sim::data_type a, sim::data_type b) {

std::cout << "Error between Initial and Final Energy" << (a-b)/b*100 << "%" << std::endl;

}
