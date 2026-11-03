#include "basic.h"
#include "constant.h"
#include "readfile.h"
#include "initialization.h"
#include "error.h"
#include "measure.h"
#include "sort_vector.h"
#include "RNA_recovery.h"

using namespace CBS_NS;

const int NBBMax = 150;
const double delETA = 20.0;
const double delDIST = 1.5;
const double delTHETA = 15.0;
const double RMSD_SMALL = 0.05;
const double RMSD_WARNING = 1.0;

RNARecovery::RNARecovery(Cobmms *cbs) : Pointers(cbs), rna_aa(cbs) {}

RNARecovery::RNARecovery(const RNARecovery &other) : Pointers(other.cbs), rna_aa(other.cbs) {
  rna_aa=other.rna_aa;
  Vbackbones_AA_atoms=other.Vbackbones_AA_atoms;
  Vbackbones_CG_xyz=other.Vbackbones_CG_xyz;
  Vbackbones_CG_N_xyz=other.Vbackbones_CG_N_xyz;
  Vetas=other.Vetas;
  Vdistance_SpPns=other.Vdistance_SpPns;
  Vdistance_SBNs=other.Vdistance_SBNs;
  Vthetas=other.Vthetas;
  Mbase_AA_atoms=other.Mbase_AA_atoms;
  Mbase_CG_xyz=other.Mbase_CG_xyz;
}

RNARecovery& RNARecovery::operator =(const RNARecovery &other) {
  if(this != &other) {
    Pointers::operator=(other.cbs);
    rna_aa=other.rna_aa;
    Vbackbones_AA_atoms=other.Vbackbones_AA_atoms;
    Vbackbones_CG_xyz=other.Vbackbones_CG_xyz;
    Vbackbones_CG_N_xyz=other.Vbackbones_CG_N_xyz;
    Vetas=other.Vetas;
    Vdistance_SpPns=other.Vdistance_SpPns;
    Vdistance_SBNs=other.Vdistance_SBNs;
    Vthetas=other.Vthetas;
    Mbase_AA_atoms=other.Mbase_AA_atoms;
    Mbase_CG_xyz=other.Mbase_CG_xyz;;
  }

  return *this;
}

RNARecovery::~RNARecovery() {
  rna_aa.clear();
  Vbackbones_AA_atoms.clear();
  Vbackbones_AA_atoms.shrink_to_fit();
  Vbackbones_CG_xyz.clear();
  Vbackbones_CG_xyz.shrink_to_fit();
  Vbackbones_CG_N_xyz.clear();
  Vbackbones_CG_N_xyz.shrink_to_fit();
  Vetas.clear();
  Vetas.shrink_to_fit();
  Vdistance_SpPns.clear();
  Vdistance_SpPns.shrink_to_fit();
  Vdistance_SBNs.clear();
  Vdistance_SBNs.shrink_to_fit();
  Vthetas.clear();
  Vthetas.shrink_to_fit();
  Mbase_AA_atoms.clear();
  Mbase_CG_xyz.clear();
}

void RNARecovery::recover_RNA_to_AA(const Molecule &rna_cg) {
  rna_aa.clear();

  RNARecovery::read_backbone_conformer();
  RNARecovery::set_base_conformer();
  
  int nthread=init->nthread;
  for(size_t ii=0;ii<rna_cg.chains.size();ii++) {
    vector<vector<Atom> > Vresidue_atoms;
    vector<vector<double> > Vresidue_backbone_xyz;
    size_t nres=rna_cg.chains[ii].residues.size();
    for(size_t jj=0;jj<nres;jj++) {
      string chaind_one=rna_cg.chains[ii].residues[jj].chainid;
      int resid_one=rna_cg.chains[ii].residues[jj].resid;
      string resname_one=rna_cg.chains[ii].residues[jj].resname;
      if(resname_one[0]=='c') resname_one=resname_one.substr(1);
      vector<double> Vxyz_ref;
      size_t natoms=rna_cg.chains[ii].residues[jj].atoms.size();
      for(size_t kk=natoms-3;kk<natoms;kk++) {
        vector<double> xyz=rna_cg.chains[ii].residues[jj].atoms[kk].get_atom_xyz();
        Vxyz_ref.insert(Vxyz_ref.end(),xyz.begin(),xyz.end());
      }
      if(Mbase_AA_atoms.count(resname_one)==0) {
        string str="RNA_recovery.cpp: Cannot find atoms for resname "+resname_one;
        error->one(str);
      }
      vector<double> Vxyz_tar=Mbase_CG_xyz[resname_one];
      vector<Atom> base_j=Mbase_AA_atoms[resname_one];
      vector<double> Vxyz_AA;
      for(size_t kk=0;kk<base_j.size();kk++) {
        base_j[kk].chainid=chaind_one;
        base_j[kk].resid=resid_one;
        base_j[kk].resname=resname_one;
        vector<double> xyz=base_j[kk].get_atom_xyz();
        Vxyz_AA.insert(Vxyz_AA.end(),xyz.begin(),xyz.end());
      }
      Superpose sp_one(Vxyz_ref,Vxyz_tar);
      sp_one.calculate_rotation_matrix();
      sp_one.reset_x2(Vxyz_AA);
      Vxyz_AA.clear();
      Vxyz_AA=sp_one.rotate_x2();
      for(size_t kk=0;kk<base_j.size();kk++) {
        base_j[kk].x=Vxyz_AA[kk*3];
        base_j[kk].y=Vxyz_AA[kk*3+1];
        base_j[kk].z=Vxyz_AA[kk*3+2];
      }
      Vresidue_atoms.push_back(base_j);

      vector<double> Vxyz_bb;
      if(jj>0) {
        Atom atom_Sp=rna_cg.chains[ii].residues[jj-1].find_atom("S");
        vector<double> xyz=atom_Sp.get_atom_xyz();
        Vxyz_bb.insert(Vxyz_bb.end(),xyz.begin(),xyz.end());
      }
      for(size_t kk=0;kk<natoms-3;kk++) {
        vector<double> xyz=rna_cg.chains[ii].residues[jj].atoms[kk].get_atom_xyz();
        Vxyz_bb.insert(Vxyz_bb.end(),xyz.begin(),xyz.end());
      }
      if(jj<nres-1) {
        Atom atom_Pn=rna_cg.chains[ii].residues[jj+1].find_atom("P");
        vector<double> xyz=atom_Pn.get_atom_xyz();
        Vxyz_bb.insert(Vxyz_bb.end(),xyz.begin(),xyz.end());
      }
      
      Atom atom_NB;
      for(size_t kk=0;kk<base_j.size();kk++) {
        if(base_j[0].resname=="A" || base_j[0].resname=="G") {
          if(base_j[kk].atomname=="N9") {
            atom_NB=base_j[kk];
            break;
          }
        }
        else if(base_j[0].resname=="C" || base_j[0].resname=="U") {
          if(base_j[kk].atomname=="N1") {
            atom_NB=base_j[kk];
            break;
          }
        }
      }
      vector<double> xyz=atom_NB.get_atom_xyz();
      Vxyz_bb.insert(Vxyz_bb.end(),xyz.begin(),xyz.end());
      
      Vresidue_backbone_xyz.push_back(Vxyz_bb);
    }

    for(size_t jj=0;jj<nres;jj++) {
      string residue_name_j=Vresidue_atoms[jj][0].chainid+"/"+to_string(Vresidue_atoms[jj][0].resid)+"/"+Vresidue_atoms[jj][0].resname;
      vector<size_t> Vcandidate_indexs;
      vector<vector<double> > VVtar_xyz;
      size_t len=Vresidue_backbone_xyz[jj].size();
      if(jj==0) {
        size_t lstart=0,lend=0;
        int flag_P=0;
        if(len==15) {
          lstart=3;
          lend=15;
          flag_P=1;
        }
        else if(len==12 && nres>1) {
          lstart=6;
          lend=15;
          flag_P=0;
        }
        else if(len==12 && nres==1) {
          lstart=3;
          lend=12;
          flag_P=1;
        }
        else if(len==9 && nres==1) {
          lstart=6;
          lend=12;
          flag_P=0;
        }
        else {
          string str="RNA_recovery.cpp: Cannot recover backbone for residue "+residue_name_j;
          str+=" --> Number_of_bead = "+to_string(len/3);
          error->one(str);
        }
        vector<double> xyz_P,xyz_S,xyz_S2,xyz_BN;
        double theta_j=180+2*delTHETA;
        double dist_j;
        if(flag_P) {
          for(size_t kk=0;kk<3;kk++)
            xyz_P.push_back(Vresidue_backbone_xyz[jj][kk]);
          for(size_t kk=3;kk<6;kk++)
            xyz_S.push_back(Vresidue_backbone_xyz[jj][kk]);
          for(size_t kk=6;kk<9;kk++)
            xyz_S2.push_back(Vresidue_backbone_xyz[jj][kk]);
          for(size_t kk=9;kk<12;kk++)
            xyz_BN.push_back(Vresidue_backbone_xyz[jj][kk]);
          theta_j=Measure::angle(xyz_P,xyz_S,xyz_S2);
          dist_j=Measure::distance(xyz_S,xyz_BN);
        }
        else {
          for(size_t kk=0;kk<3;kk++)
            xyz_S.push_back(Vresidue_backbone_xyz[jj][kk]);
          for(size_t kk=6;kk<9;kk++)
            xyz_BN.push_back(Vresidue_backbone_xyz[jj][kk]);
          dist_j=Measure::distance(xyz_S,xyz_BN);
        }

        vector<double> Vdiff;
        if(flag_P) {
          for(size_t kk=0;kk<Vdistance_SBNs.size();kk++) {
            double ddist=fabs(Vdistance_SBNs[kk]-dist_j);
            double dtheta_1=fabs(Vthetas[kk]-theta_j);
            double dtheta_2=fabs(Vthetas[kk]+180-theta_j);
            double dtheta_3=fabs(Vthetas[kk]-180-theta_j);
            double dtheta=min(dtheta_1,dtheta_2);
            dtheta=min(dtheta,dtheta_3);
            double diff=ddist/delDIST+dtheta/delTHETA;
            Vdiff.push_back(diff);
          }
        }
        else {
          for(size_t kk=0;kk<Vdistance_SBNs.size();kk++) {
            double ddist=fabs(Vdistance_SBNs[kk]-dist_j);
            double diff=ddist/delDIST;
            Vdiff.push_back(diff);
          } 
        }
        vector<size_t> Vsort=SortVector::sort_vector_ascending(Vdiff);
        int nconf=(int)Vsort.size();
        if(nconf>NBBMax) nconf=NBBMax;
        for(int kk=0;kk<nconf;kk++) {
          size_t index_k=Vsort[kk];
          vector<double> Vtar_xyz;
          for(size_t ll=lstart;ll<lend;ll++)
            Vtar_xyz.push_back(Vbackbones_CG_xyz[index_k][ll]);
          for(size_t ll=0;ll<3;ll++)
            Vtar_xyz.push_back(Vbackbones_CG_N_xyz[index_k][ll]);
          VVtar_xyz.push_back(Vtar_xyz);
          Vcandidate_indexs.push_back(index_k);
        }
      }
      else if(nres>1 && jj==nres-1) {
        vector<double> xyz_P,xyz_S,xyz_S2,xyz_BN;
        for(size_t kk=3;kk<6;kk++)
          xyz_P.push_back(Vresidue_backbone_xyz[jj][kk]);
        for(size_t kk=6;kk<9;kk++)
          xyz_S.push_back(Vresidue_backbone_xyz[jj][kk]);
        for(size_t kk=9;kk<12;kk++)
          xyz_S2.push_back(Vresidue_backbone_xyz[jj][kk]);
        for(size_t kk=12;kk<15;kk++)
          xyz_BN.push_back(Vresidue_backbone_xyz[jj][kk]);
        double theta_j=Measure::angle(xyz_P,xyz_S,xyz_S2);
        double dist_j=Measure::distance(xyz_S,xyz_BN);
        
        vector<double> Vdiff;
        for(size_t kk=0;kk<Vdistance_SBNs.size();kk++) {
          double ddist=fabs(Vdistance_SBNs[kk]-dist_j);
          double dtheta_1=fabs(Vthetas[kk]-theta_j);
          double dtheta_2=fabs(Vthetas[kk]+180-theta_j);
          double dtheta_3=fabs(Vthetas[kk]-180-theta_j);
          double dtheta=min(dtheta_1,dtheta_2);
          dtheta=min(dtheta,dtheta_3);
          double diff=ddist/delDIST+dtheta/delTHETA;
          Vdiff.push_back(diff);
        }
        vector<size_t> Vsort=SortVector::sort_vector_ascending(Vdiff);
        int nconf=(int)Vsort.size();
        if(nconf>NBBMax) nconf=NBBMax;
        for(int kk=0;kk<nconf;kk++) {
          size_t index_k=Vsort[kk];
          vector<double> Vtar_xyz;
          for(size_t ll=0;ll<12;ll++)
            Vtar_xyz.push_back(Vbackbones_CG_xyz[index_k][ll]);
          for(size_t ll=0;ll<3;ll++)
            Vtar_xyz.push_back(Vbackbones_CG_N_xyz[index_k][ll]);
          VVtar_xyz.push_back(Vtar_xyz);
          Vcandidate_indexs.push_back(index_k);
        }
      }
      else {
        vector<double> Vdiff;
        vector<double> xyz_Sp,xyz_P,xyz_S,xyz_S2,xyz_Pn,xyz_BN;
        for(size_t kk=0;kk<3;kk++)
          xyz_Sp.push_back(Vresidue_backbone_xyz[jj][kk]);
        for(size_t kk=3;kk<6;kk++)
          xyz_P.push_back(Vresidue_backbone_xyz[jj][kk]);
        for(size_t kk=6;kk<9;kk++)
          xyz_S.push_back(Vresidue_backbone_xyz[jj][kk]);
        for(size_t kk=9;kk<12;kk++)
          xyz_S2.push_back(Vresidue_backbone_xyz[jj][kk]);
        for(size_t kk=12;kk<15;kk++)
          xyz_Pn.push_back(Vresidue_backbone_xyz[jj][kk]);
        for(size_t kk=15;kk<18;kk++)
          xyz_BN.push_back(Vresidue_backbone_xyz[jj][kk]);
        double eta_j=Measure::dihedral(xyz_Sp,xyz_P,xyz_S,xyz_Pn);
        double dist1_j=Measure::distance(xyz_Sp,xyz_Pn);
        double dist2_j=Measure::distance(xyz_S,xyz_BN);
        double theta_j=Measure::angle(xyz_P,xyz_S,xyz_S2);
        for(size_t kk=0;kk<Vetas.size();kk++) {
          double deta_1=fabs(Vetas[kk]-eta_j);
          double deta_2=fabs(Vetas[kk]+360-eta_j);
          double deta_3=fabs(Vetas[kk]-360-eta_j);
          double deta=min(deta_1,deta_2);
          deta=min(deta,deta_3);
          if(deta>delETA) continue;
          double ddist1=fabs(Vdistance_SpPns[kk]-dist1_j);
          if(ddist1>delDIST) continue;
          double ddist2=fabs(Vdistance_SBNs[kk]-dist2_j);
          if(ddist2>delDIST) continue;
          double dtheta_1=fabs(Vthetas[kk]-theta_j);
          double dtheta_2=fabs(Vthetas[kk]+180-theta_j);
          double dtheta_3=fabs(Vthetas[kk]-180-theta_j);
          double dtheta=min(dtheta_1,dtheta_2);
          dtheta=min(dtheta,dtheta_3);
          if(dtheta>delTHETA) continue;
          double diff=ddist1/delDIST+ddist2/delDIST+deta/delETA+dtheta/delTHETA;
          Vdiff.push_back(diff);
        }
        if(Vdiff.size()==0) {
          for(size_t kk=0;kk<Vetas.size();kk++) {
            double deta_1=fabs(Vetas[kk]-eta_j);
            double deta_2=fabs(Vetas[kk]+360-eta_j);
            double deta_3=fabs(Vetas[kk]-360-eta_j);
            double deta=min(deta_1,deta_2);
            deta=min(deta,deta_3);
            double ddist1=fabs(Vdistance_SpPns[kk]-dist1_j);
            double ddist2=fabs(Vdistance_SBNs[kk]-dist2_j);
            double dtheta_1=fabs(Vthetas[kk]-theta_j);
            double dtheta_2=fabs(Vthetas[kk]+180-theta_j);
            double dtheta_3=fabs(Vthetas[kk]-180-theta_j);
            double dtheta=min(dtheta_1,dtheta_2);
            dtheta=min(dtheta,dtheta_3);
            double diff=ddist1/delDIST+ddist2/delDIST+deta/delETA+dtheta/delTHETA;
            Vdiff.push_back(diff);
          }
        }
        vector<size_t> Vsort=SortVector::sort_vector_ascending(Vdiff);
        int nconf=(int)Vsort.size();
        if(nconf>NBBMax) nconf=NBBMax;
        for(int kk=0;kk<nconf;kk++) {
          size_t index_k=Vsort[kk];
          vector<double> Vtar_xyz=Vbackbones_CG_xyz[index_k];
          for(size_t ll=0;ll<3;ll++)
            Vtar_xyz.push_back(Vbackbones_CG_N_xyz[index_k][ll]);
          VVtar_xyz.push_back(Vtar_xyz);
          Vcandidate_indexs.push_back(index_k);
        }
      }
        
      vector<Superpose> Vsp;
      for(int nn=0;nn<nthread;nn++) {
        Superpose sp_one;
        sp_one.set_x1(Vresidue_backbone_xyz[jj]);
        Vsp.push_back(sp_one);
      }
      vector<double> Vrmsds;
      size_t index_min;
      double rmsd_min=MathConst::BIG;
      int flag_break=0;
      int ncal=(int)VVtar_xyz.size();
      int nloop=ncal/nthread;
      int nlast=ncal%nthread;
      for(int nl=0;nl<nloop;nl++) {
        int nl0=nl*nthread;
        vector<vector<double> > Vxyz_loop;
        Vxyz_loop.resize(nthread);
        for(int ni=0;ni<nthread;ni++)
          Vxyz_loop[ni]=VVtar_xyz[nl0+ni];
        vector<double> rmsd_loop;
        rmsd_loop.resize(nthread);
        vector<std::thread> threads;
        for(int ni=0;ni<nthread;ni++)
          threads.emplace_back(std::thread(RNARecovery::calculate_RMSD,std::ref(Vsp[ni]),std::ref(Vxyz_loop[ni]),std::ref(rmsd_loop),ni));
        for(int ni=0;ni<nthread;ni++)
          threads[ni].join();
        for(int ni=0;ni<nthread;ni++) {
          Vrmsds.push_back(rmsd_loop[ni]);
          if(rmsd_loop[ni]<=RMSD_SMALL) {
            rmsd_min=rmsd_loop[ni];
            index_min=Vrmsds.size()-1;
            flag_break=1;
            break;
          }
        }
        if(flag_break) break;
      }
      if(nlast!=0 && flag_break==0) {
        int nl0=nloop*nthread;
        vector<vector<double> > Vxyz_last;
        Vxyz_last.resize(nlast);
        for(int ni=0;ni<nlast;ni++)
          Vxyz_last[ni]=VVtar_xyz[nl0+ni];
        vector<double> rmsd_last;
        rmsd_last.resize(nlast);
        vector<std::thread> threads;
        for(int ni=0;ni<nlast;ni++)
          threads.emplace_back(std::thread(RNARecovery::calculate_RMSD,std::ref(Vsp[ni]),std::ref(Vxyz_last[ni]),std::ref(rmsd_last),ni));
        for(int ni=0;ni<nlast;ni++)
          threads[ni].join();
        for(int ni=0;ni<nlast;ni++) {
          Vrmsds.push_back(rmsd_last[ni]);
          if(rmsd_last[ni]<=RMSD_SMALL) {
            rmsd_min=rmsd_last[ni];
            index_min=Vrmsds.size()-1;
            flag_break=1;
            break;
          }
        }
      }

      if(flag_break==0) {
        for(size_t kk=0;kk<Vrmsds.size();kk++) {
          if(Vrmsds[kk]<rmsd_min) {
            rmsd_min=Vrmsds[kk];
            index_min=kk;
          }
        }
      }
        
      Superpose sp_one;
      sp_one.set_x1(Vresidue_backbone_xyz[jj]);
      sp_one.set_x2(VVtar_xyz[index_min]);
      rmsd_min=sp_one.get_RMSD();
      vector<Atom> backbone_j=Vbackbones_AA_atoms[Vcandidate_indexs[index_min]];
      vector<double> Vtar_AA;
      for(size_t kk=0;kk<backbone_j.size();kk++) {
        backbone_j[kk].chainid=Vresidue_atoms[jj][0].chainid;
        backbone_j[kk].resid=Vresidue_atoms[jj][0].resid;
        backbone_j[kk].resname=Vresidue_atoms[jj][0].resname;
        vector<double> xyz=backbone_j[kk].get_atom_xyz();
        Vtar_AA.insert(Vtar_AA.end(),xyz.begin(),xyz.end());
      }
      Vtar_AA=sp_one.rotate_x2(Vtar_AA);
      for(size_t kk=0;kk<backbone_j.size();kk++) {
        backbone_j[kk].x=Vtar_AA[3*kk];
        backbone_j[kk].y=Vtar_AA[3*kk+1];
        backbone_j[kk].z=Vtar_AA[3*kk+2];
      }
      vector<Atom> residue_j=backbone_j;
      residue_j.insert(residue_j.end(),Vresidue_atoms[jj].begin(),Vresidue_atoms[jj].end());
      Vresidue_atoms[jj]=residue_j;
      if(rmsd_min>RMSD_WARNING) {
        string str="RNA_recovery.cpp: RMSD= "+to_string(rmsd_min)+" for residue "+residue_name_j;
        error->warning(str);
      }
    }

    Chain chain_one(cbs);
    chain_one.chain_label=Chain::rna;
    for(size_t jj=0;jj<nres;jj++) {
      Residue residue_one(cbs);
      residue_one.set_residue(Vresidue_atoms[jj]);
      chain_one.insert_residue(residue_one,Chain::end);
    }
    rna_aa.insert_chain(chain_one);
  }

  int index_one=0;
  for(size_t ii=0;ii<rna_aa.chains.size();ii++)
    index_one=rna_aa.chains[ii].set_atomid(index_one);
}

void RNARecovery::recover_RNA_to_AA_v2(const Molecule &rna_cg) {
  rna_aa.clear();

  RNARecovery::read_backbone_conformer();
  RNARecovery::set_base_conformer();

  for(size_t ii=0;ii<Vbackbones_CG_xyz.size();ii++) {
    vector<double> Txyz=Vbackbones_CG_xyz[ii];
    vector<double> Vxyz;
    for(size_t jj=0;jj<Txyz.size();jj++) {
      if(jj>=9 && jj<=11) continue;  //remove xyz_S2
      Vxyz.push_back(Txyz[jj]);
    }
    Vbackbones_CG_xyz[ii]=Vxyz;
  }
  
  int nthread=init->nthread;
  for(size_t ii=0;ii<rna_cg.chains.size();ii++) {
    vector<vector<Atom> > Vresidue_atoms;
    vector<vector<double> > Vresidue_backbone_xyz;
    size_t nres=rna_cg.chains[ii].residues.size();
    for(size_t jj=0;jj<nres;jj++) {
      string chaind_one=rna_cg.chains[ii].residues[jj].chainid;
      int resid_one=rna_cg.chains[ii].residues[jj].resid;
      string resname_one=rna_cg.chains[ii].residues[jj].resname;
      if(resname_one[0]=='c') resname_one=resname_one.substr(1);
      vector<double> Vxyz_ref;
      size_t natoms=rna_cg.chains[ii].residues[jj].atoms.size();
      for(size_t kk=natoms-3;kk<natoms;kk++) {
        vector<double> xyz=rna_cg.chains[ii].residues[jj].atoms[kk].get_atom_xyz();
        Vxyz_ref.insert(Vxyz_ref.end(),xyz.begin(),xyz.end());
      }
      if(Mbase_AA_atoms.count(resname_one)==0) {
        string str="RNA_recovery.cpp: Cannot find atoms for resname "+resname_one;
        error->one(str);
      }
      vector<double> Vxyz_tar=Mbase_CG_xyz[resname_one];
      vector<Atom> base_j=Mbase_AA_atoms[resname_one];
      vector<double> Vxyz_AA;
      for(size_t kk=0;kk<base_j.size();kk++) {
        base_j[kk].chainid=chaind_one;
        base_j[kk].resid=resid_one;
        base_j[kk].resname=resname_one;
        vector<double> xyz=base_j[kk].get_atom_xyz();
        Vxyz_AA.insert(Vxyz_AA.end(),xyz.begin(),xyz.end());
      }
      Superpose sp_one(Vxyz_ref,Vxyz_tar);
      sp_one.calculate_rotation_matrix();
      sp_one.reset_x2(Vxyz_AA);
      Vxyz_AA.clear();
      Vxyz_AA=sp_one.rotate_x2();
      for(size_t kk=0;kk<base_j.size();kk++) {
        base_j[kk].x=Vxyz_AA[kk*3];
        base_j[kk].y=Vxyz_AA[kk*3+1];
        base_j[kk].z=Vxyz_AA[kk*3+2];
      }
      Vresidue_atoms.push_back(base_j);

      vector<double> Vxyz_bb;
      if(jj>0) {
        Atom atom_Sp=rna_cg.chains[ii].residues[jj-1].find_atom("S");
        vector<double> xyz=atom_Sp.get_atom_xyz();
        Vxyz_bb.insert(Vxyz_bb.end(),xyz.begin(),xyz.end());
      }
      for(size_t kk=0;kk<natoms-3;kk++) {
        vector<double> xyz=rna_cg.chains[ii].residues[jj].atoms[kk].get_atom_xyz();
        Vxyz_bb.insert(Vxyz_bb.end(),xyz.begin(),xyz.end());
      }
      if(jj<nres-1) {
        Atom atom_Pn=rna_cg.chains[ii].residues[jj+1].find_atom("P");
        vector<double> xyz=atom_Pn.get_atom_xyz();
        Vxyz_bb.insert(Vxyz_bb.end(),xyz.begin(),xyz.end());
      }
    
      Atom atom_NB;
      for(size_t kk=0;kk<base_j.size();kk++) {
        if(base_j[0].resname=="A" || base_j[0].resname=="G") {
          if(base_j[kk].atomname=="N9") {
            atom_NB=base_j[kk];
            break;
          }
        }
        else if(base_j[0].resname=="C" || base_j[0].resname=="U") {
          if(base_j[kk].atomname=="N1") {
            atom_NB=base_j[kk];
            break;
          }
        }
      }
      vector<double> xyz=atom_NB.get_atom_xyz();
      Vxyz_bb.insert(Vxyz_bb.end(),xyz.begin(),xyz.end());
      
      Vresidue_backbone_xyz.push_back(Vxyz_bb);
    }

    for(size_t jj=0;jj<nres;jj++) {
      string residue_name_j=Vresidue_atoms[jj][0].chainid+"/"+to_string(Vresidue_atoms[jj][0].resid)+"/"+Vresidue_atoms[jj][0].resname;
      vector<size_t> Vcandidate_indexs;
      vector<vector<double> > VVtar_xyz;
      size_t len=Vresidue_backbone_xyz[jj].size();
      if(jj==0) {
        size_t lstart=0,lend=0;
        int flag_P=0;
        if(len==12) {
          lstart=3;
          lend=12;
          flag_P=1;
        }
        else if(len==9 && nres>1) {
          lstart=6;
          lend=12;
          flag_P=0;
        }
        else if(len==9 && nres==1) {
          lstart=3;
          lend=9;
          flag_P=1;
        }
        //else if(len==6 && nres==1) {
        //  lstart=6;
        //  lend=9;
        //  flag_P=0;
        //}
        else {
          string str="RNA_recovery.cpp: Cannot recover backbone for residue "+residue_name_j;
          str+=" --> Number_of_bead = "+to_string(len/3);
          error->one(str);
        }
        vector<double> xyz_P,xyz_S,xyz_BN;
        double dist_j;
        if(flag_P) {
          for(size_t kk=0;kk<3;kk++)
            xyz_P.push_back(Vresidue_backbone_xyz[jj][kk]);
          for(size_t kk=3;kk<6;kk++)
            xyz_S.push_back(Vresidue_backbone_xyz[jj][kk]);
          for(size_t kk=6;kk<9;kk++)
            xyz_BN.push_back(Vresidue_backbone_xyz[jj][kk]);
          dist_j=Measure::distance(xyz_S,xyz_BN);
        }
        else {
          for(size_t kk=0;kk<3;kk++)
            xyz_S.push_back(Vresidue_backbone_xyz[jj][kk]);
          for(size_t kk=3;kk<6;kk++)
            xyz_BN.push_back(Vresidue_backbone_xyz[jj][kk]);
          dist_j=Measure::distance(xyz_S,xyz_BN);
        }

        vector<double> Vdiff;
        if(flag_P) {
          for(size_t kk=0;kk<Vdistance_SBNs.size();kk++) {
            double ddist=fabs(Vdistance_SBNs[kk]-dist_j);
            double diff=ddist/delDIST;
            Vdiff.push_back(diff);
          }
        }
        else {
          for(size_t kk=0;kk<Vdistance_SBNs.size();kk++) {
            double ddist=fabs(Vdistance_SBNs[kk]-dist_j);
            double diff=ddist/delDIST;
            Vdiff.push_back(diff);
          } 
        }
        vector<size_t> Vsort=SortVector::sort_vector_ascending(Vdiff);
        int nconf=(int)Vsort.size();
        if(nconf>NBBMax) nconf=NBBMax;
        for(int kk=0;kk<nconf;kk++) {
          size_t index_k=Vsort[kk];
          vector<double> Vtar_xyz;
          for(size_t ll=lstart;ll<lend;ll++)
            Vtar_xyz.push_back(Vbackbones_CG_xyz[index_k][ll]);
          for(size_t ll=0;ll<3;ll++)
            Vtar_xyz.push_back(Vbackbones_CG_N_xyz[index_k][ll]);
          VVtar_xyz.push_back(Vtar_xyz);
          Vcandidate_indexs.push_back(index_k);
        }
      }
      else if(nres>1 && jj==nres-1) {
        vector<double> xyz_P,xyz_S,xyz_BN;
        for(size_t kk=3;kk<6;kk++)
          xyz_P.push_back(Vresidue_backbone_xyz[jj][kk]);
        for(size_t kk=6;kk<9;kk++)
          xyz_S.push_back(Vresidue_backbone_xyz[jj][kk]);
        for(size_t kk=9;kk<12;kk++)
          xyz_BN.push_back(Vresidue_backbone_xyz[jj][kk]);
        double dist_j=Measure::distance(xyz_S,xyz_BN);
        
        vector<double> Vdiff;
        for(size_t kk=0;kk<Vdistance_SBNs.size();kk++) {
          double ddist=fabs(Vdistance_SBNs[kk]-dist_j);
          double diff=ddist/delDIST;
          Vdiff.push_back(diff);
        }
        vector<size_t> Vsort=SortVector::sort_vector_ascending(Vdiff);
        int nconf=(int)Vsort.size();
        if(nconf>NBBMax) nconf=NBBMax;
        for(int kk=0;kk<nconf;kk++) {
          size_t index_k=Vsort[kk];
          vector<double> Vtar_xyz;
          for(size_t ll=0;ll<9;ll++)
            Vtar_xyz.push_back(Vbackbones_CG_xyz[index_k][ll]);
          for(size_t ll=0;ll<3;ll++)
            Vtar_xyz.push_back(Vbackbones_CG_N_xyz[index_k][ll]);
          VVtar_xyz.push_back(Vtar_xyz);
          Vcandidate_indexs.push_back(index_k);
        }
      }
      else {
        vector<double> Vdiff;
        vector<double> xyz_Sp,xyz_P,xyz_S,xyz_Pn,xyz_BN;
        for(size_t kk=0;kk<3;kk++)
          xyz_Sp.push_back(Vresidue_backbone_xyz[jj][kk]);
        for(size_t kk=3;kk<6;kk++)
          xyz_P.push_back(Vresidue_backbone_xyz[jj][kk]);
        for(size_t kk=6;kk<9;kk++)
          xyz_S.push_back(Vresidue_backbone_xyz[jj][kk]);
        for(size_t kk=9;kk<12;kk++)
          xyz_Pn.push_back(Vresidue_backbone_xyz[jj][kk]);
        for(size_t kk=12;kk<15;kk++)
          xyz_BN.push_back(Vresidue_backbone_xyz[jj][kk]);
        double eta_j=Measure::dihedral(xyz_Sp,xyz_P,xyz_S,xyz_Pn);
        double dist1_j=Measure::distance(xyz_Sp,xyz_Pn);
        double dist2_j=Measure::distance(xyz_S,xyz_BN);
        for(size_t kk=0;kk<Vetas.size();kk++) {
          double deta_1=fabs(Vetas[kk]-eta_j);
          double deta_2=fabs(Vetas[kk]+360-eta_j);
          double deta_3=fabs(Vetas[kk]-360-eta_j);
          double deta=min(deta_1,deta_2);
          deta=min(deta,deta_3);
          if(deta>delETA) continue;
          double ddist1=fabs(Vdistance_SpPns[kk]-dist1_j);
          if(ddist1>delDIST) continue;
          double ddist2=fabs(Vdistance_SBNs[kk]-dist2_j);
          if(ddist2>delDIST) continue;
          double diff=ddist1/delDIST+ddist2/delDIST+deta/delETA;
          Vdiff.push_back(diff);
        }
        if(Vdiff.size()==0) {
          for(size_t kk=0;kk<Vetas.size();kk++) {
            double deta_1=fabs(Vetas[kk]-eta_j);
            double deta_2=fabs(Vetas[kk]+360-eta_j);
            double deta_3=fabs(Vetas[kk]-360-eta_j);
            double deta=min(deta_1,deta_2);
            deta=min(deta,deta_3);
            double ddist1=fabs(Vdistance_SpPns[kk]-dist1_j);
            double ddist2=fabs(Vdistance_SBNs[kk]-dist2_j);
            double diff=ddist1/delDIST+ddist2/delDIST+deta/delETA;
            Vdiff.push_back(diff);
          }
        }
        vector<size_t> Vsort=SortVector::sort_vector_ascending(Vdiff);
        int nconf=(int)Vsort.size();
        if(nconf>NBBMax) nconf=NBBMax;
        for(int kk=0;kk<nconf;kk++) {
          size_t index_k=Vsort[kk];
          vector<double> Vtar_xyz=Vbackbones_CG_xyz[index_k];
          for(size_t ll=0;ll<3;ll++)
            Vtar_xyz.push_back(Vbackbones_CG_N_xyz[index_k][ll]);
          VVtar_xyz.push_back(Vtar_xyz);
          Vcandidate_indexs.push_back(index_k);
        }
      }
        
      vector<Superpose> Vsp;
      for(int nn=0;nn<nthread;nn++) {
        Superpose sp_one;
        sp_one.set_x1(Vresidue_backbone_xyz[jj]);
        Vsp.push_back(sp_one);
      }
      vector<double> Vrmsds;
      size_t index_min;
      double rmsd_min=MathConst::BIG;
      int flag_break=0;
      int ncal=(int)VVtar_xyz.size();
      int nloop=ncal/nthread;
      int nlast=ncal%nthread;
      for(int nl=0;nl<nloop;nl++) {
        int nl0=nl*nthread;
        vector<vector<double> > Vxyz_loop;
        Vxyz_loop.resize(nthread);
        for(int ni=0;ni<nthread;ni++)
          Vxyz_loop[ni]=VVtar_xyz[nl0+ni];
        vector<double> rmsd_loop;
        rmsd_loop.resize(nthread);
        vector<std::thread> threads;
        for(int ni=0;ni<nthread;ni++)
          threads.emplace_back(std::thread(RNARecovery::calculate_RMSD,std::ref(Vsp[ni]),std::ref(Vxyz_loop[ni]),std::ref(rmsd_loop),ni));
        for(int ni=0;ni<nthread;ni++)
          threads[ni].join();
        for(int ni=0;ni<nthread;ni++) {
          Vrmsds.push_back(rmsd_loop[ni]);
          if(rmsd_loop[ni]<=RMSD_SMALL) {
            rmsd_min=rmsd_loop[ni];
            index_min=Vrmsds.size()-1;
            flag_break=1;
            break;
          }
        }
        if(flag_break) break;
      }
      if(nlast!=0 && flag_break==0) {
        int nl0=nloop*nthread;
        vector<vector<double> > Vxyz_last;
        Vxyz_last.resize(nlast);
        for(int ni=0;ni<nlast;ni++)
          Vxyz_last[ni]=VVtar_xyz[nl0+ni];
        vector<double> rmsd_last;
        rmsd_last.resize(nlast);
        vector<std::thread> threads;
        for(int ni=0;ni<nlast;ni++)
          threads.emplace_back(std::thread(RNARecovery::calculate_RMSD,std::ref(Vsp[ni]),std::ref(Vxyz_last[ni]),std::ref(rmsd_last),ni));
        for(int ni=0;ni<nlast;ni++)
          threads[ni].join();
        for(int ni=0;ni<nlast;ni++) {
          Vrmsds.push_back(rmsd_last[ni]);
          if(rmsd_last[ni]<=RMSD_SMALL) {
            rmsd_min=rmsd_last[ni];
            index_min=Vrmsds.size()-1;
            flag_break=1;
            break;
          }
        }
      }

      if(flag_break==0) {
        for(size_t kk=0;kk<Vrmsds.size();kk++) {
          if(Vrmsds[kk]<rmsd_min) {
            rmsd_min=Vrmsds[kk];
            index_min=kk;
          }
        }
      }
        
      Superpose sp_one;
      sp_one.set_x1(Vresidue_backbone_xyz[jj]);
      sp_one.set_x2(VVtar_xyz[index_min]);
      rmsd_min=sp_one.get_RMSD();
      vector<Atom> backbone_j=Vbackbones_AA_atoms[Vcandidate_indexs[index_min]];
      vector<double> Vtar_AA;
      for(size_t kk=0;kk<backbone_j.size();kk++) {
        backbone_j[kk].chainid=Vresidue_atoms[jj][0].chainid;
        backbone_j[kk].resid=Vresidue_atoms[jj][0].resid;
        backbone_j[kk].resname=Vresidue_atoms[jj][0].resname;
        vector<double> xyz=backbone_j[kk].get_atom_xyz();
        Vtar_AA.insert(Vtar_AA.end(),xyz.begin(),xyz.end());
      }
      Vtar_AA=sp_one.rotate_x2(Vtar_AA);
      for(size_t kk=0;kk<backbone_j.size();kk++) {
        backbone_j[kk].x=Vtar_AA[3*kk];
        backbone_j[kk].y=Vtar_AA[3*kk+1];
        backbone_j[kk].z=Vtar_AA[3*kk+2];
      }
      vector<Atom> residue_j=backbone_j;
      residue_j.insert(residue_j.end(),Vresidue_atoms[jj].begin(),Vresidue_atoms[jj].end());
      Vresidue_atoms[jj]=residue_j;
      if(rmsd_min>RMSD_WARNING) {
        string str="RNA_recovery.cpp: RMSD= "+to_string(rmsd_min)+" for residue "+residue_name_j;
        error->warning(str);
      }
    }

    Chain chain_one(cbs);
    chain_one.chain_label=Chain::rna;
    for(size_t jj=0;jj<nres;jj++) {
      Residue residue_one(cbs);
      residue_one.set_residue(Vresidue_atoms[jj]);
      chain_one.insert_residue(residue_one,Chain::end);
    }
    rna_aa.insert_chain(chain_one);
  }

  int index_one=0;
  for(size_t ii=0;ii<rna_aa.chains.size();ii++)
    index_one=rna_aa.chains[ii].set_atomid(index_one);
}

void RNARecovery::clear() {
  rna_aa.clear();
  Vbackbones_AA_atoms.clear();
  Vbackbones_AA_atoms.shrink_to_fit();
  Vbackbones_CG_xyz.clear();
  Vbackbones_CG_xyz.shrink_to_fit();
  Vbackbones_CG_N_xyz.clear();
  Vbackbones_CG_N_xyz.shrink_to_fit();
  Vetas.clear();
  Vetas.shrink_to_fit();
  Vdistance_SpPns.clear();
  Vdistance_SpPns.shrink_to_fit();
  Vdistance_SBNs.clear();
  Vdistance_SBNs.shrink_to_fit();
  Vthetas.clear();
  Vthetas.shrink_to_fit();
  Mbase_AA_atoms.clear();
  Mbase_CG_xyz.clear();
}

void RNARecovery::read_backbone_conformer() {
  Vbackbones_AA_atoms.clear();
  Vbackbones_CG_xyz.clear();
  Vbackbones_CG_N_xyz.clear();
  Vetas.clear();
  Vdistance_SpPns.clear();
  Vdistance_SBNs.clear();
  Vthetas.clear();

  size_t Nbackbone_atoms=12;

  string fstr=init->INPUT_DIR+"RNA_backbone.dat";
  vector<string> Vlines=read_file(fstr);
  vector<Atom> backbone;
  for(size_t ii=0;ii<Vlines.size();ii++) {
    string &sline=Vlines[ii];
    if(sline.size()<2) continue;
    if(sline.find("#")!=string::npos && sline.find("Eta")!=string::npos) {
      vector<string> stoken=split_string(sline,":");
      if(stoken.size()!=8) {
        string str="RNA_recovery.cpp: Wrong head line in read RNA_backbone.dat --> "+sline;
        error->one(str);
      }
      Vetas.push_back(stod(stoken[2]));
      Vdistance_SpPns.push_back(stod(stoken[4]));
      Vdistance_SBNs.push_back(stod(stoken[5]));
      Vthetas.push_back(stod(stoken[7]));
      if(backbone.size()>=Nbackbone_atoms) {
        vector<Atom> Tatoms;
        for(size_t kk=0;kk<Nbackbone_atoms;kk++)
          Tatoms.push_back(backbone[kk]);
        Vbackbones_AA_atoms.push_back(Tatoms);
        vector<double> Vxyz,Vxyz_N;
        for(size_t kk=Nbackbone_atoms;kk<backbone.size();kk++) {
          vector<double> xyz=backbone[kk].get_atom_xyz();
          if(backbone[kk].atomname!="BN") 
            Vxyz.insert(Vxyz.end(),xyz.begin(),xyz.end());
          else 
            Vxyz_N.insert(Vxyz_N.end(),xyz.begin(),xyz.end());
        }
        Vbackbones_CG_xyz.push_back(Vxyz);
        Vbackbones_CG_N_xyz.push_back(Vxyz_N);
      }
      backbone.clear();
      continue;
    }
    vector<string> stoken=split_string(sline);
    if(stoken.size()!=4) {
      string str="RNA_recovery.cpp: Wrong line in read RNA_backbone.dat --> "+sline;
      error->one(str);
    }
    string atomname_one=stoken[0];
    double x_one=stod(stoken[1]);
    double y_one=stod(stoken[2]);
    double z_one=stod(stoken[3]);
    Atom atom_one;
    atom_one.set_atom(atomname_one,x_one,y_one,z_one);
    backbone.push_back(atom_one);
  }
  if(backbone.size()>=Nbackbone_atoms) {
    vector<Atom> Tatoms;
    for(size_t kk=0;kk<Nbackbone_atoms;kk++)
      Tatoms.push_back(backbone[kk]);
    Vbackbones_AA_atoms.push_back(Tatoms);
    vector<double> Vxyz,Vxyz_N;
    for(size_t kk=Nbackbone_atoms;kk<backbone.size();kk++) {
      vector<double> xyz=backbone[kk].get_atom_xyz();
      if(backbone[kk].atomname!="BN")
        Vxyz.insert(Vxyz.end(),xyz.begin(),xyz.end());
      else
        Vxyz_N.insert(Vxyz_N.end(),xyz.begin(),xyz.end());
    }
    Vbackbones_CG_xyz.push_back(Vxyz);
    Vbackbones_CG_N_xyz.push_back(Vxyz_N);
    backbone.clear();
  }

  for(size_t ii=0;ii<Vbackbones_AA_atoms.size();ii++) {
    if(Vbackbones_AA_atoms[ii].size()!=Nbackbone_atoms) {
      string str="RNA_recovery.cpp: Number of atoms in backbone conformer is not 12 --> ";
      str+=to_string(Vbackbones_AA_atoms[ii].size());
      error->one(str);
    }
  }

  if(Vbackbones_AA_atoms.size()!=Vetas.size()) {
    string str="RNA_recovery.cpp: Number of backbone conformers is not matched --> ";
    str+=to_string(Vbackbones_AA_atoms.size())+" / "+to_string(Vetas.size());
    error->one(str);
  }
}

void RNARecovery::set_base_conformer() {
  Mbase_AA_atoms.clear();
  Mbase_CG_xyz.clear();

  vector<Atom> &baseA=init->RESIDUES_XYZ["baseA"];
  Mbase_AA_atoms.emplace("A",baseA);
  vector<Atom> &baseG=init->RESIDUES_XYZ["baseG"];
  Mbase_AA_atoms.emplace("G",baseG);
  vector<Atom> &baseC=init->RESIDUES_XYZ["baseC"];
  Mbase_AA_atoms.emplace("C",baseC);
  vector<Atom> &baseU=init->RESIDUES_XYZ["baseU"];
  Mbase_AA_atoms.emplace("U",baseU);

  vector<Atom> &cA=init->RESIDUES_XYZ["cA"];
  vector<double> Vxyz;
  for(size_t ii=3;ii<cA.size();ii++) {
    vector<double> xyz=cA[ii].get_atom_xyz();
    Vxyz.insert(Vxyz.end(),xyz.begin(),xyz.end());
  }
  Mbase_CG_xyz.emplace("A",Vxyz);
  vector<Atom> &cG=init->RESIDUES_XYZ["cG"];
  Vxyz.clear();
  for(size_t ii=3;ii<cG.size();ii++) {
    vector<double> xyz=cG[ii].get_atom_xyz();
    Vxyz.insert(Vxyz.end(),xyz.begin(),xyz.end());
  }
  Mbase_CG_xyz.emplace("G",Vxyz);
  vector<Atom> &cC=init->RESIDUES_XYZ["cC"];
  Vxyz.clear();
  for(size_t ii=3;ii<cC.size();ii++) {
    vector<double> xyz=cC[ii].get_atom_xyz();
    Vxyz.insert(Vxyz.end(),xyz.begin(),xyz.end());
  }
  Mbase_CG_xyz.emplace("C",Vxyz);
  vector<Atom> &cU=init->RESIDUES_XYZ["cU"];
  Vxyz.clear();
  for(size_t ii=3;ii<cU.size();ii++) {
    vector<double> xyz=cU[ii].get_atom_xyz();
    Vxyz.insert(Vxyz.end(),xyz.begin(),xyz.end());
  }
  Mbase_CG_xyz.emplace("U",Vxyz);
}

int RNARecovery::calculate_RMSD(Superpose &sp_one, const vector<double> &Vtar, vector<double> &Vrmsd, int nth) {
  sp_one.set_x2(Vtar);
  Vrmsd[nth]=sp_one.get_RMSD();

  return 0;
}