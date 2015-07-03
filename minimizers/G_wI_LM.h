#ifndef XFEL_LEGACY_QUADRANTS_H
#define XFEL_LEGACY_QUADRANTS_H
#include <scitbx/array_family/flex_types.h>
#include <scitbx/array_family/shared.h>
#include <map>
#include <scitbx/lstbx/normal_equations.h>
//#include <ctime>
//#include <iostream>

namespace xscale {
  class xscale_G_wI_LM: public scitbx::lstbx::normal_equations::non_linear_ls<double> {

    public:
      xscale_G_wI_LM(int n_parameters):
        scitbx::lstbx::normal_equations::non_linear_ls<double>(n_parameters),
        jacobian_one_row(n_parameters)
        {}
      void set_cpp_data(int x) {
       residuals = scitbx::af::shared<double>(x);
      }

      inline
      void
      build_up_normal_matrix(bool objective_only,
                                    const shared_double& I,
                                    const shared_double& w,
                                    const shared_size_t& hkl,
                                    const shared_size_t& frames,
                                    const shared_double& S,
                                    const shared_double& G,
                                    const shared_double& Ih,
                                    const int n_frames) {
        double Gm=0;
        double f=0;
        int m=0;
        int h=0;
        double dtdIh=0;
        double t = 0;
        double residuals=0;
        using namespace std;
//        clock_t start, end;
//        float avg=0;
//        int counter=0;
        for (int i=0;i<I.size();i++) {
//          counter+=1;
          m=frames[i];
          h = hkl[i];
          residuals = w[i]*(I[i] - (G[m] * Ih[h] * S[h]));
          jacobian_one_row[m] = -w[i] * Ih[h] * S[h];
          jacobian_one_row[h+n_frames] = -w[i] * G[m] * S[h];
//          start = clock();
          add_equation(residuals, jacobian_one_row.const_ref(), 1.);
//          end = clock();
//          avg+=(float)((end - start) / (float)CLOCKS_PER_SEC);
          jacobian_one_row[m] = 0;
          jacobian_one_row[h+n_frames]= 0;
        }
//      cout<<"The average time for add equation is: "<<avg/counter<<" sec"<<endl;
    }

/*      inline
      void
      build_up_sparse_normal_matrix(bool objective_only,
                                    const shared_double& I,
                                    const shared_double& w,
                                    const shared_size_t& hkl,
                                    const shared_size_t& frames,
                                    const shared_double& S,
                                    const shared_double& G,
                                    const shared_double& Ih,
                                    const int n_frames) {

        double Gm=0;
        double f=0;
        int m=0;
        int h=0;
        double dtdIh=0;
        double t = 0;
        double residuals=0;
//        using namespace std;
//        clock_t start, end;
//        float avg=0;
//        int counter=0;
        for (int i=0;i<I.size();i++) {
//          counter+=1;
          std::map<int, double> row;
          m=frames[i];
          h = hkl[i];
          residuals = w[i]*(I[i] - (G[m] * Ih[h] * S[h]));
          row[m] = -w[i] * Ih[h] * S[h];
          row[h+n_frames] = -w[i] * G[m] * S[h];
//          start = clock();
          add_sparse_equation(residuals, row);
//          end = clock();
//          avg+=(float)((end - start) / (float)CLOCKS_PER_SEC);
      }
//      cout<<"The average time for add equation is: "<<avg/counter<<" sec"<<endl;
     }

*/
    private:
      scitbx::af::shared<double> jacobian_one_row, residuals;
  };
}

#endif
