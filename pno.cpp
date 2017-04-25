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
    }

    // commit the result to energy
    this->set_value(energy, energy_);
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











// //**********************************************************************
// // Implementation of MP2
// //**********************************************************************


// // This is a basic implementation of (iterative) MP2 energy.
// // This is formulated in terms of occupied and unoccupied states represented
// // as LCAOs, hence it's derived from LCAOWfn, aka ::mpqc::lcao::LCAOWavefunction
// // .
// // This class can only compute the energy, this is indicated by deriving it from
// // Provides<Energy> (this introduces virtual methods can_evaluate() and
// // evaluate() that specify, respectively, the ability to compute the energy and
// // how the energy is computed).
// class MP2 : public LCAOWfn, public Provides<Energy> {
//  public:
//   // a few abbreviations to make the code less verbose
//   using Array = TA::DistArray<TA::TensorD, TA::SparsePolicy>;
//   using RHF = lcao::RHF<TA::TensorD, TA::SparsePolicy>;

//   // the KeyVal constructor takes a KeyVal object that represents a keyword
//   // group of the input file that corresponds to this object (see the "mp2"
//   // group in the mp2.json file that accompanies this example).
//   // The Keyval object will be queried for all keywords needed by
//   // the KeyVal ctor of LCAOWfn, as well as keyword "ref" that specifies
//   // the reference wave function.
//   MP2(const KeyVal& kv) : LCAOWfn(kv) {
//     ref_wfn_ = kv.class_ptr<RHF>("ref");
//     if (!ref_wfn_)
//       throw InputError("missing reference RHF wave function", __FILE__,
//                        __LINE__, "ref");
//   }

//  private:
//   // This implements the Energy::Provider::can_evaluate() virtual
//   // function. This function returns true is the energy object can be computed.
//   bool can_evaluate(Energy* energy) override {
//     // can only compute energies (not forces (energy->order() == 1),
//     // hessians, or higher derivatives)
//     return energy->order() == 0;
//   }

//   // This implements the Energy::Provider::evaluate() virtual function.
//   // This function computes the MP2 energy and assigns it to the Energy object.
//   void evaluate(Energy* energy) override {
//     // how precisely to compute the energy
//     auto target_precision = energy->target_precision(0);

//     // if has not been computed to the desired precision
//     // (or never computed at all) ...
//     if (computed_precision_ > target_precision) {
//       // compute reference to higher precision than this wfn
//       auto target_ref_precision = target_precision / 100.;
//       auto ref_energy =
//           std::make_shared<Energy>(ref_wfn_, target_ref_precision);
//       ref_wfn_->evaluate(ref_energy.get());

//       // use the reference orbitals to populate the orbital space registry
//       init_sdref(ref_wfn_, target_ref_precision);

//       // this actually computes the energy
//       double mp2_corr_energy = compute_mp2_energy(target_precision);

//       energy_ = ref_energy->energy() + mp2_corr_energy;
//     }

//     // commit the result to energy
//     this->set_value(energy, energy_);
//   }

//   // This function actually solves the MP1 equations and returns the MP2 energy.
//   double compute_mp2_energy(double target_precision) {
//     // fac is an LCAOFactory object which evaluates integrals in terms of AOs
//     // and LCAOs
//     auto& fac = this->lcao_factory();
//     auto& world = fac.world();
//     // ofac is an OrbitalSpaceRegistry that defines the LCAO spaces that fac can
//     // use ofac was populated by the init_sdref() call above
//     auto& ofac = fac.orbital_registry();

//     auto nocc = ofac.retrieve("m").rank();
//     ExEnv::out0() << "nocc = " << nocc;
//     auto nocc_act = ofac.retrieve("i").rank();
//     auto nvir = ofac.retrieve("a").rank();
//     auto nfzc = nocc - nocc_act;

//     auto F = fac.compute(L"(p|F|q)");
//     Eigen::VectorXd eps_p = array_ops::array_to_eigen(F).diagonal();
//     // replicated diagonal elements of Fo
//     auto eps_o = eps_p.segment(nfzc, nocc_act);
//     for (int i=0; i<eps_o.size(); ++i) {
//       std::cout << "eps(i,i) = " << eps_o[i] << std::endl;
//     }
//     // replicated diagonal elements of Fv
//     auto eps_v = eps_p.tail(nvir);

//     // G_iajb
//     auto G = fac.compute(L"(i a|G|j b)");
//     // Fij
//     auto Fo = fac.compute(L"(i|F|j)");
//     // Fab
//     auto Fv = fac.compute(L"(a|F|b)");

//     //K_aibj
//     auto K = fac.compute(L"(a i|G|b j)");
//     const auto ktrange = K.trange();

//     std::vector<Eigen::MatrixXd> K_vec(nocc_act*nocc_act); 

//       for (int i=0; i<nocc_act; ++i) {
//         for (int j=0; j<nocc_act; ++j) {
//           std::array<int,4> idx = {{0,i,0,j}};
//           auto tile_ij = ktrange.element_to_tile(idx);
//           auto ord_ij = ktrange.tiles_range().ordinal(tile_ij);
//           TA::TensorD K_ij = K.find(ord_ij);
//           auto ext = K_ij.range().extent_data();
//           Eigen::MatrixXd K_ij_mat = TA::eigen_map(K_ij, ext[0]*ext[1], ext[2]*ext[3]);
//           K_vec[i*nocc_act + j] = K_ij_mat;
//           //ExEnv::out0() << "Hello, World!" << std::endl;
//           //ExEnv::out0() << "T^" << i << "," << j << ":\n" << T_ij_mat << std::endl;
//         }
//       }


//       std::vector<Eigen::MatrixXd> T_vec(nocc_act*nocc_act);

//       for (int i=0; i<nocc_act; ++i) {
//         auto eps_i = eps_o[i];
//         for (int j=0; j<nocc_act; ++j) {
//           auto eps_j = eps_o[j];
//           Eigen::MatrixXd K_ij = K_vec[i*nocc_act + j];
//           Eigen::MatrixXd T_ij(nvir, nvir);
//           for (int a=0; a<nvir; ++a) {
//             auto eps_a = eps_v[a];
//             for (int b=0; b<nvir; ++b) {
//               auto eps_b = eps_v[b];
//               T_ij(a,b) = -K_ij(a,b)/(eps_a + eps_b - eps_i - eps_j);
//               //T_vec[i*nocc_act + j] = T_ij;
//             }
//           }
//           T_vec[i*nocc_act + j] = T_ij;
//         }
//       }


//       std::vector<Eigen::MatrixXd> T_tilde_vec(nocc_act*nocc_act);

//       for (int i=0; i<nocc_act; ++i) {
//         for (int j=0; j<nocc_act; ++j) {
//           Eigen::MatrixXd T_ij = T_vec[i*nocc_act + j];
//           Eigen::MatrixXd T_ji = T_vec[j*nocc_act + i];
//           Eigen::MatrixXd T_tilde_ij = 2*T_ij - T_ji;
//           T_tilde_vec[i*nocc_act + j] = T_tilde_ij;
//         }
//       }

//       std::vector<Eigen::MatrixXd> D_vec(nocc_act*nocc_act);

//       for (int i=0; i<nocc_act; ++i) {
//         for (int j=0; j<nocc_act; ++j) {
//           auto delta_ij = (i == j)? 1 : 0;
//           Eigen::MatrixXd T_ij = T_vec[i*nocc_act + j];
//           Eigen::MatrixXd T_tilde_ij = T_tilde_vec[i*nocc_act + j];
//           Eigen::MatrixXd D_ij = (T_tilde_ij.adjoint()*T_ij + T_tilde_ij*T_ij.adjoint())
//           /(1.0 + delta_ij);
//           D_vec[i*nocc_act + j] = D_ij;
//         }
//       }




//     // zero out amplitudes
//     if (!T_.is_initialized()) {
//       T_ = Array(world, G.trange(), G.shape());
//       T_.fill(0.0);
//     }

//     // lambda function will be used to do a Jacobi update of the residual
//     auto jacobi_update = [eps_o, eps_v](TA::TensorD& result_tile) {

//       const auto& range = result_tile.range();
//       double norm = 0.0;
//       for (const auto& i : range) {
//         const auto result_abij = result_tile[i] / (eps_o[i[0]] - eps_v[i[1]] +
//                                                    eps_o[i[2]] - eps_v[i[3]]);
//         result_tile[i] = result_abij;
//         norm += result_abij * result_abij;
//       }
//       return std::sqrt(norm);
//     };

//     // solve the MP1 equations
//     auto converged = false;
//     auto iter = 0;
//     auto energy = +1.0;
//     ExEnv::out0() << "Start solving MP2 Energy\n" << std::endl;
//     while (not converged) {
//       Array R;
//       R("i,a,j,b") = G("i,a,j,b") + Fv("a,c") * T_("i,c,j,b") +
//                      Fv("b,c") * T_("i,a,j,c") - Fo("i,k") * T_("k,a,j,b") -
//                      Fo("j,k") * T_("i,a,k,b");

//       // const auto trange = G.trange();
//       // std::array<int,4> idx{{3,0,2,0}};
//       // auto tile_with_elem_3121 = trange.element_to_tile(idx);
//       // auto tord = trange.tiles_range().ordinal(tile_with_elem_3121);
//       // ExEnv::out0() << "Tile " << tile_with_elem_3121 << " contains elem 3,0,2,0.\n";
//       // ExEnv::out0() << "Tile " << tord << " contains elem 3,0,2,0.\n";
//       // TA::TensorD tile = G.find(tord);
//       // auto ext = tile.range().extent_data();
//       // Eigen::MatrixXd tile_mat = TA::eigen_map(tile,ext[0] * ext[1], ext[2] * ext[3]);
//       //ExEnv::out0() << "i=3 j=2_{ab} = \n" << tile_mat << "\n";

//       // estimate the MP2 energy ... the Hylleraas formula is quadratic in error
//       double updated_energy =
//           (G("i,a,j,b") + R("i,a,j,b")).dot(2 * T_("i,a,j,b") - T_("i,b,j,a"));
//       ExEnv::out0() << indent << "Iteration: " << iter
//                     << " Energy: " << updated_energy << std::endl;

//       computed_precision_ = std::abs(updated_energy - energy);
//       energy = updated_energy;

//       // update the amplitudes, if needed
//       converged = computed_precision_ <= target_precision;
//       if (not converged) {
//         // R^{ij}_{ab} -> R^{ij}_{ab} / (F^i_i + F^j_j - F^a_a - F^b_b)
//         TA::foreach_inplace(R, jacobi_update);
//         // need a fence here since foreach_inplace mutates the contents of R
//         // as a side effect.
//         // N.B. most TiledArray ops will not need a fence (except to control
//         //      the resource use)
//         world.gop.fence();
//         T_("i,a,j,b") += R("i,a,j,b");
//         ++iter;
//       }
//     }

//     return energy;
//   }



//   // This reimplements the ::mpqc::Wavefunction::obsolete() virtual function.
//   // It gets called when, for example, the atomic positions get updated in
//   // geometry optimization.
//   void obsolete() override {
//     LCAOWfn::obsolete();
//     ref_wfn_->obsolete();
//     computed_precision_ = std::numeric_limits<double>::max();
//     energy_ = 0.0;
//   }

//   std::shared_ptr<RHF> ref_wfn_;
//   Array T_;
//   double computed_precision_ = std::numeric_limits<double>::max();
//   double energy_ = 0.0;
// };

// // This macro registers the KeyVal constructor of our MP2 class and associates
// // it with the "MP2" key, so that the KeyVal class knows how to find it when it
// // finds this key as the object type in the user input.
// MPQC_CLASS_EXPORT2("MP2", MP2);

// // Creating this variable forces the code for the MP2 class to be linked into
// // the mp2 executable (otherwise the MPQC main function will not see any
// // references to this class and thus the linker will simply skip it).
// mpqc::detail::ForceLink<MP2> fl;



// //**********************************************************************
// // Implementation of Pair Natural Orbitals (PNOs)
// //**********************************************************************

// // class PNOs : public LCAOWfn, public Provides<Energy> {
// //   public:
// //   // a few abbreviations to make the code less verbose
// //   using Array = TA::DistArray<TA::TensorD, TA::SparsePolicy>;
// //   using RHF = lcao::RHF<TA::TensorD, TA::SparsePolicy>;

//   // Check to make sure that reference RHF wavefunction exists
//   // MP2(const KeyVal& kv) : LCAOWfn(kv) {
//   //   ref_wfn_ = kv.class_ptr<RHF>("ref");
//   //   if (!ref_wfn_)
//   //     throw InputError("missing reference RHF wave function", __FILE__,
//   //                      __LINE__, "ref");
//   // }

//   // Declare 4-d array to hold T amplitudes

//   // Construct tile boundary vectors for each of the 4 dimensions

//   // Dimensions 0 and 1 are of size nocc and holds nocc tiles, each with size 1
//   // std::vector<std::size_t> occ_tile_boundaries;
//   // for (std::size_t i=0; i<= nocc; i += 1) {
//   //   occ_tile_boundaries.push_back(i);
//   // }

//   // // Dimensions 2 and 3 are of size nvirt and each holds hold a single tile
//   // std::vector<std::size_t> virt_tile_boundaries;
//   // for (std::size_t i=0; i<= nvirt; i += nvirt) {
//   //   virt_tile_boundaries.push_back(i);
//   // }
  
//   // // Construct a set of 1D TiledRanges
//   // std::vector<>
//   // // Each tile corresponds to a matrix for 4-index integrals with a single
//   // // i and j value but all possible a and b values
//   // T = Array(world, )



// // };



//     //Eigen::MatrixXd fock = TA::array_to_eigen(F);
//     //std::cout << "fock:\n" << fock << std::endl; 
//     //TA::Array<double, 4> ao_ints;


//     // // Construct tile boundary vector for dimension 0 (occ)
//     // std::vector<int> tile_boundaries0(nocc+1);
//     // for (int i=0; i<= nocc; ++i) {
//     //   //tile_boundaries0.push_back(i);
//     //   tile_boundaries0[i] = i;
//     // }

//     // // Construct tile boundary vector for dimension 1 (occ)
//     // std::vector<int> tile_boundaries1(nocc+1);
//     // for (int i=0; i<=nocc; ++i) {
//     //   //tile_boundaries1.push_back(i);
//     //   tile_boundaries1[i] = 1;
//     // }

//     // // Construct tile boundary vector for dimension 2 (unocc)
//     // std::vector<int> tile_boundaries2(1);
//     // for (int i=0; i<=nvir; ++nvir) {
//     //   //tile_boundaries2.push_back(i);
//     //   tile_boundaries2[i] = i;
//     // }

//     // // Construct tile boundary vector for dimension 3 (unocc)
//     // std::vector<int> tile_boundaries3(1);
//     // for (int i=0; i<=nvir; ++nvir) {
//     //   //tile_boundaries3.push_back(i);
//     //   tile_boundaries3[i] = i;
//     // }

//     // for (int i=0; i<=tile_boundaries0.size(); ++i) {
//     //   std::cout << tile_boundaries0[i] << std::endl;
//     // }

//     // Construct a set of 1D TiledRanges
//     // std::vector<TA::TiledRange1> ranges(4, (TA::TiledRange1(tile_boundaries0.begin(),
//     // tile_boundaries0.end()), TA::TiledRange1(tile_boundaries1.begin(),
//     // tile_boundaries1.end()), TA::TiledRange1(tile_boundaries2.begin(),
//     // tile_boundaries2.end()), TA::TiledRange1(tile_boundaries3.begin(),
//     // tile_boundaries3.end())));

//     // // Construct the 4D TiledRange
//     // TA::TiledRange trange(ranges.begin(), ranges.end());

//     // //std::cout << "trange = " << trange << std::endl;

//     // TA::Array<double, 4> ao_ints(world, trange);





//   // double compute_pnos(double target_precision) {
//   //   // fac is an LCAOFactory object which evaluates integrals in terms of AOs
//   //   // and LCAOs
//   //   auto& fac = this->lcao_factory();
//   //   auto& world = fac.world();
//   //   // ofac is an OrbitalSpaceRegistry that defines the LCAO spaces that fac can
//   //   // use ofac was populated by the init_sdref() call above
//   //   auto& ofac = fac.orbital_registry();

//   //   auto nocc = ofac.retrieve("m").rank();
//   //   ExEnv::out0() << "nocc = " << nocc;
//   //   auto nocc_act = ofac.retrieve("i").rank();
//   //   auto nvir = ofac.retrieve("a").rank();
//   //   auto nfzc = nocc - nocc_act;

//   //   auto F = fac.compute(L"(p|F|q)");
//   //   Eigen::VectorXd eps_p = array_ops::array_to_eigen(F).diagonal();
//   //   // replicated diagonal elements of Fo
//   //   auto eps_o = eps_p.segment(nfzc, nocc_act);
//   //   // replicated diagonal elements of Fv
//   //   auto eps_v = eps_p.tail(nvir);

//   //   // G_iajb
//   //   //auto G = fac.compute(L"(i a|G|j b)");

//   //   // K_aibj
//   //   auto K = fac.compute(L"(a i|G|b j)");
//   //   const auto ktrange = K.trange();

//   //   std::vector<Eigen::MatrixXd> T_vec(nocc_act*nocc_act); 

//   //    for (int i=0; i<nocc_act; ++i) {
//   //      for (int j=0; j<nocc_act; ++j) {
//   //        std::array<int,4> idx = {{0,i,0,j}};
//   //        auto tile_ij = ktrange.element_to_tile(idx);
//   //        auto ord_ij = ktrange.tiles_range().ordinal(tile_ij);
//   //        TA::TensorD T_ij = K.find(ord_ij);
//   //        auto ext = T_ij.range().extent_data();
//   //        Eigen::MatrixXd T_ij_mat = TA::eigen_map(T_ij, ext[0]*ext[1], ext[2]*ext[3]);
//   //        T_vec[i*nocc_act + j] = T_ij_mat;
//   //        ExEnv::out0() << "Hello, World!" << std::endl;
//   //        //ExEnv::out0() << "T^" << i << "," << j << ":\n" << T_ij_mat << std::endl;
//   //      }
//   //    }



