#include "mpqc/chemistry/qc/lcao/scf/rhf.h"
#include "mpqc/chemistry/qc/lcao/wfn/lcao_wfn.h"
#include "mpqc/chemistry/qc/properties/energy.h"
#include "mpqc/util/keyval/forcelink.h"

using namespace mpqc;

using LCAOWfn = lcao::LCAOWavefunction<TA::TensorD, TA::SparsePolicy>;


//**********************************************************************
// Implementation of PNOs using MP2 amplitudes
//**********************************************************************

class PNO : public LCAOWfn, public Provides<Energy> {
  public:
    // A few abbreviations to make the code less verbose
    using Array = TA::DistArray<TA::TensorD, TA::SparsePolicy>;
    using RHF = lcao::RHF<TA::TensorD, TA::SparsePolicy>;



  // the KeyVal constructor takes a KeyVal object that represents a keyword
  // group of the input file that corresponds to this object (see the "mp2"
  // group in the mp2.json file that accompanies this example).
  // The Keyval object will be queried for all keywords needed by
  // the KeyVal ctor of LCAOWfn, as well as keyword "ref" that specifies
  // the reference wave function.
  PNO(const KeyVal& kv) : LCAOWfn(kv) {
    ref_wfn_ = kv.class_ptr<RHF>("ref");
    if (!ref_wfn_)
      throw InputError("missing reference RHF wave function", __FILE__,
                       __LINE__, "ref");
  }

  
  private: 
    bool can_evaluate(Energy* energy) override {
    // can only compute energies (not forces (energy->order() == 1),
    // hessians, or higher derivatives)
    return energy->order() == 0;
    }

    // This implements the Energy::Provider::evaluate() virtual function.
    // This function computes the MP2 energy and assigns it to the Energy object.
    void evaluate(Energy* energy) override {
      // how precisely to compute the energy
      auto target_precision = energy->target_precision(0);

      // if has not been computed to the desired precision
      // (or never computed at all) ...
      if (computed_precision_ > target_precision) {
        // compute reference to higher precision than this wfn
        auto target_ref_precision = target_precision / 100.;
        auto ref_energy =
          std::make_shared<Energy>(ref_wfn_, target_ref_precision);
        ref_wfn_->evaluate(ref_energy.get());

        // use the reference orbitals to populate the orbital space registry
        init_sdref(ref_wfn_, target_ref_precision);

        // this actually computes the energy
        double mp2_corr_energy = compute_mp2_energy(target_precision);

        energy_ = ref_energy->energy() + mp2_corr_energy;

        double threshold = 1.0e-7;

        int pno = compute_pnos(target_precision, threshold);
    }

    // commit the result to energy
    this->set_value(energy, energy_);
  }

    int compute_pnos(double target_precision, double threshold) {
    //std::vector<Eigen::MatrixXd> make_D(double target_precision) {

      auto& fac = this->lcao_factory();
      auto& world = fac.world();
      auto& ofac = fac.orbital_registry();

      auto nocc = ofac.retrieve("m").rank();
      ExEnv::out0() << "nocc = " << nocc;
      auto nocc_act = ofac.retrieve("i").rank();
      auto nvir = ofac.retrieve("a").rank();
      auto nfzc = nocc - nocc_act;

      auto F = fac.compute(L"(p|F|q)");
      Eigen::VectorXd eps_p = array_ops::array_to_eigen(F).diagonal();
      auto eps_o = eps_p.segment(nfzc, nocc_act);
      auto eps_v = eps_p.tail(nvir);

      // K_aibj
      auto K = fac.compute(L"(a i|G|b j)");
      const auto ktrange = K.trange();

      std::vector<Eigen::MatrixXd> K_vec(nocc_act*nocc_act); 

      for (int i=0; i<nocc_act; ++i) {
        for (int j=0; j<nocc_act; ++j) {
          std::array<int,4> tile_ij = {{0,i,0,j}};
          auto ord_ij = ktrange.tiles_range().ordinal(tile_ij);
          ExEnv::out0() << "[" << tile_ij[0] << "," << tile_ij[1]
          << "," << tile_ij[2] << "," << tile_ij[3] << "] occurs at " << ord_ij <<std::endl;
          TA::TensorD K_ij = K.find(ord_ij);
          auto ext = K_ij.range().extent_data();
          Eigen::MatrixXd K_ij_mat = TA::eigen_map(K_ij, ext[0]*ext[1], ext[2]*ext[3]);
          K_vec[i*nocc_act + j] = K_ij_mat;
        }
      }

      std::vector<Eigen::MatrixXd> T_vec(nocc_act*nocc_act);

      for (int i=0; i<nocc_act; ++i) {
        auto eps_i = eps_o[i];
        for (int j=0; j<nocc_act; ++j) {
          auto eps_j = eps_o[j];
          Eigen::MatrixXd K_ij = K_vec[i*nocc_act + j];
          Eigen::MatrixXd T_ij(nvir, nvir);
          for (int a=0; a<nvir; ++a) {
            auto eps_a = eps_v[a];
            for (int b=0; b<nvir; ++b) {
              auto eps_b = eps_v[b];
              T_ij(a,b) = -K_ij(a,b)/(eps_a + eps_b - eps_i - eps_j);
              //T_vec[i*nocc_act + j] = T_ij;
            }
          }
          T_vec[i*nocc_act + j] = T_ij;
        }
      }

      std::vector<Eigen::MatrixXd> T_tilde_vec(nocc_act*nocc_act);

      for (int i=0; i<nocc_act; ++i) {
        for (int j=0; j<nocc_act; ++j) {
          Eigen::MatrixXd T_ij = T_vec[i*nocc_act + j];
          Eigen::MatrixXd T_ji = T_vec[j*nocc_act + i];
          Eigen::MatrixXd T_tilde_ij = 2*T_ij - T_ji;
          T_tilde_vec[i*nocc_act + j] = T_tilde_ij;
        }
      }

      std::vector<Eigen::MatrixXd> D_vec(nocc_act*nocc_act);

      for (int i=0; i<nocc_act; ++i) {
        for (int j=0; j<nocc_act; ++j) {
          auto delta_ij = (i == j)? 1 : 0;
          Eigen::MatrixXd T_ij = T_vec[i*nocc_act + j];
          Eigen::MatrixXd T_tilde_ij = T_tilde_vec[i*nocc_act + j];
          Eigen::MatrixXd D_ij = (T_tilde_ij.adjoint()*T_ij + T_tilde_ij*T_ij.adjoint())
          /(1.0 + delta_ij);
          D_vec[i*nocc_act + j] = D_ij;
        }
      }
  
      std::vector<Eigen::MatrixXd> pno_vec(nocc_act*nocc_act);
      std::vector<Eigen::VectorXd> occ_vec(nocc_act*nocc_act);
      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;

      for (int i=0; i<nocc_act; ++i) {
        for (int j=0; j<nocc_act; ++j) {
          Eigen::MatrixXd D_ij = D_vec[i*nocc_act + j];
          es.compute(D_ij);
          Eigen::VectorXd occ = es.eigenvalues();
          Eigen::MatrixXd pnos = es.eigenvectors();
          Eigen::VectorXd occ_keep(occ.size());
          Eigen::MatrixXd pnos_keep(pnos.rows(),pnos.cols());
          //Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, 0, nocc_act, 1> occ_keep;
          //Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, 0, nocc_act, nocc_act> pnos_keep;
          int idx=0;
          for (int k=0; k<occ.size(); ++k) {
            if (occ[k] > threshold) {
              occ_keep[idx] = occ[k];
              pnos_keep.col(idx) = pnos.col(k);
              idx += 1;
            }
            else
              continue;
          }
          occ_vec[i*nocc_act + j] = occ_keep;
          pno_vec[i*nocc_act + j] = pnos_keep;
          ExEnv::out0() << "For (" << i << "," << j << ") there are " << idx
          << " PNOs with occupation numbers greater than " << threshold << std::endl;
        }
      }
      return 0;
    } 



    double compute_mp2_energy(double target_precision) {
      double energy = 0.05;
      return energy;
    }


  // This reimplements the ::mpqc::Wavefunction::obsolete() virtual function.
  // It gets called when, for example, the atomic positions get updated in
  // geometry optimization.
  void obsolete() override {
    LCAOWfn::obsolete();
    ref_wfn_->obsolete();
    computed_precision_ = std::numeric_limits<double>::max();
    energy_ = 0.0;
  }

  std::shared_ptr<RHF> ref_wfn_;
  //Array T_;
  double computed_precision_ = std::numeric_limits<double>::max();
  double energy_ = 0.0;


};


// This macro registers the KeyVal constructor of our MP2 class and associates
// it with the "MP2" key, so that the KeyVal class knows how to find it when it
// finds this key as the object type in the user input.
MPQC_CLASS_EXPORT2("PNO", PNO);

// Creating this variable forces the code for the MP2 class to be linked into
// the mp2 executable (otherwise the MPQC main function will not see any
// references to this class and thus the linker will simply skip it).
mpqc::detail::ForceLink<PNO> fl;




