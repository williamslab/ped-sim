// ped-sim: pedigree simulation tool
//
// This program is distributed under the terms of the GNU General Public License

#ifndef COINTERFERE_H
#define COINTERFERE_H

using namespace std;

class COInterfere {
  public:
    COInterfere(double _nu[2], double _p[2], double _length[2]) {
      for(int i = 0; i < 2; i++) {
	nu[i] = _nu[i];
	p[i] = _p[i];
	length[i] = _length[i];
      }
      initStartProb();
    }

    void simStahl(vector<double> &locations, int sex, mt19937 &randomGen);

  private:
    void initStartProb();

    static const int N_BINS4START = 10000;

    // interference parameters, index 0 for male, 1 for female
    double nu[2];
    double p[2];
    double length[2];

    // for simulating: probability of first crossover occuring at a position
    double startProb[2][N_BINS4START];
};

#endif // COINTERFERE_H
