#include "basic.h"
#include "constant.h"
#include "initialization.h"
#include "error.h"
#include "superpose.h"
#include "math_extra.h"
#include "measure.h"
#include "assemble_RNA_motif.h"

using namespace CBS_NS;
using namespace MathExtra;


const double ATOM_DISTANCE_MIN = 2.0;
const int N_CLASH_ATOMS = 20;

AssembleRNAMotif::AssembleRNAMotif(Cobmms *cbs) : Pointers(cbs), motif_adjust(cbs) {
  moved_distance = 0.0; 
}

AssembleRNAMotif::AssembleRNAMotif(const AssembleRNAMotif &other) : Pointers(other.cbs), motif_adjust(other.cbs) {
  motif_adjust=other.motif_adjust;
  overlapped_resenames=other.overlapped_resenames;
}

AssembleRNAMotif& AssembleRNAMotif::operator=(const AssembleRNAMotif &other) {
  if(this != &other) {
    Pointers::operator=(other.cbs);
    motif_adjust=other.motif_adjust;
    overlapped_resenames=other.overlapped_resenames;
  }

  return *this;
}

AssembleRNAMotif::~AssembleRNAMotif() {
  motif_adjust.clear();
  overlapped_resenames.clear();
  overlapped_resenames.shrink_to_fit();
  moved_distance=0.0;
}

void AssembleRNAMotif::assemble_RNA_motif(const vector<RNAMotif> &Vmotif_reference, const RNAMotif &motif_target, int flag_move) {
  vector<Atom> Vatoms_ref;
  for(size_t ii=0;ii<Vmotif_reference.size();ii++) {
    vector<Atom> Vatoms_list=Vmotif_reference[ii].get_atom_list();
    Vatoms_ref.insert(Vatoms_ref.end(),Vatoms_list.begin(),Vatoms_list.end());
  }

  motif_adjust=motif_target;

  overlapped_resenames.clear();

  /*
  int flag_template=0;
  for(size_t ii=0;ii<motif_adjust.ends.size();ii++) { //if have user_templates, using user_template only as reference
    string end_name_i=motif_adjust.ends[ii];
    for(size_t jj=0;jj<Vmotif_reference.size();jj++) {
      int check=0;
      for(size_t kk=0;kk<Vmotif_reference[jj].Vresnames.size();kk++) {
        if(Vmotif_reference[jj].type.find("user_template")==string::npos) continue;
        flag_template=1;
        if(Vmotif_reference[jj].Vresnames[kk]==end_name_i) {
          overlapped_resenames.push_back(end_name_i);
          check=1;
          break;
        }
      }
      if(check) break;
    }
  }

  if(flag_template && overlapped_resenames.size()<2) {
    for(size_t ii=0;ii<motif_adjust.Vresnames.size();ii++) {
      string res_name_i=motif_adjust.Vresnames[ii];
      for(size_t jj=0;jj<Vmotif_reference.size();jj++) {
        if(Vmotif_reference[jj].type.find("user_template")==string::npos) continue;
        int check=0;
        for(size_t kk=0;kk<Vmotif_reference[jj].Vresnames.size();kk++) {
          if(Vmotif_reference[jj].Vresnames[kk]==res_name_i) {
            if(count(overlapped_resenames.begin(),overlapped_resenames.end(),res_name_i)==0)
              overlapped_resenames.push_back(res_name_i);
            check=1;
            break;
          }
        }
        if(check) break;
      }
    }
  }

  if(overlapped_resenames.size()<2) {
    for(size_t ii=0;ii<motif_adjust.ends.size();ii++) {
      string end_name_i=motif_adjust.ends[ii];
      if(count(overlapped_resenames.begin(),overlapped_resenames.end(),end_name_i)) continue;
      for(size_t jj=0;jj<Vmotif_reference.size();jj++) {
        int check=0;
        for(size_t kk=0;kk<Vmotif_reference[jj].Vresnames.size();kk++) {
          if(Vmotif_reference[jj].Vresnames[kk]==end_name_i) {
            overlapped_resenames.push_back(end_name_i);
            check=1;
            break;
          }
        }
        if(check) break;
      }
    }
  }
  */
  
  int flag_template=0;
  for(size_t ii=0;ii<Vmotif_reference.size();ii++) {
    if(Vmotif_reference[ii].type.find("user_template")==string::npos) continue; //if have user_templates, using user_template only as reference
    flag_template=1;
    for(size_t jj=0;jj<Vmotif_reference[ii].ends.size();jj++) {
      string end_name_j=Vmotif_reference[ii].ends[jj];
      for(size_t kk=0;kk<motif_adjust.Vresnames.size();kk++) {
        if(motif_adjust.Vresnames[kk]==end_name_j) {
          overlapped_resenames.push_back(end_name_j);
          break;
        }
      }
    }
  }

  if(flag_template && overlapped_resenames.size()<2) {
    for(size_t ii=0;ii<Vmotif_reference.size();ii++) {
      if(Vmotif_reference[ii].type.find("user_template")==string::npos) continue;
      for(size_t jj=0;jj<Vmotif_reference[ii].Vresnames.size();jj++) {
        string res_name_j=Vmotif_reference[ii].Vresnames[jj];
        for(size_t kk=0;kk<motif_adjust.Vresnames.size();kk++) {
          if(motif_adjust.Vresnames[kk]==res_name_j) {
            if(count(overlapped_resenames.begin(),overlapped_resenames.end(),res_name_j)==0)
              overlapped_resenames.push_back(res_name_j);
            break;
          }
        }
      }
    }
  }

  if(overlapped_resenames.size()<2) {
    for(size_t ii=0;ii<Vmotif_reference.size();ii++) {
      for(size_t jj=0;jj<Vmotif_reference[ii].ends.size();jj++) {
        string end_name_j=Vmotif_reference[ii].ends[jj];
        if(count(overlapped_resenames.begin(),overlapped_resenames.end(),end_name_j)) continue;
        for(size_t kk=0;kk<motif_adjust.Vresnames.size();kk++) {
          if(motif_adjust.Vresnames[kk]==end_name_j) {
            overlapped_resenames.push_back(end_name_j);
            break;
          }
        }
      }
    }
  }
 
  if(overlapped_resenames.size()<2) {
    for(size_t ii=0;ii<Vmotif_reference.size();ii++) {
      for(size_t jj=0;jj<Vmotif_reference[ii].Vresnames.size();jj++) {
        string res_name_j=Vmotif_reference[ii].Vresnames[jj];
        if(count(overlapped_resenames.begin(),overlapped_resenames.end(),res_name_j)) continue;
        for(size_t kk=0;kk<motif_adjust.Vresnames.size();kk++) {
          if(motif_adjust.Vresnames[kk]==res_name_j) {
            overlapped_resenames.push_back(res_name_j);
            break;
          }
        }
      }
    }
  }

  
  if(overlapped_resenames.size()==0) {
    string str="assemble_RNA_motif.cpp: Cannot find overlapped residues! --> "+motif_target.get_RNA_motif_info();
    error->one(str);
  }

  vector<double> Vxyz_ref,Vxyz_tar;
  vector<double> Vxyz_interface_tar;
  for(size_t ii=0;ii<overlapped_resenames.size();ii++) {
    int check=0;
    for(size_t jj=0;jj<Vmotif_reference.size();jj++) {
      for(size_t kk=0;kk<Vmotif_reference[jj].chains.size();kk++) {
        for(size_t ll=0;ll<Vmotif_reference[jj].chains[kk].residues.size();ll++) {
          if(Vmotif_reference[jj].chains[kk].residues[ll].get_residue_name()==overlapped_resenames[ii]) {
            vector<double> xyz=Vmotif_reference[jj].chains[kk].residues[ll].get_residue_xyz();
            Vxyz_ref.insert(Vxyz_ref.end(),xyz.begin(),xyz.end());
            check=1;
            break;
          }
        }
        if(check) break;
      }
      if(check) break;
    }
    
    check=0;
    for(size_t jj=0;jj<motif_adjust.chains.size();jj++) {
      for(size_t kk=0;kk<motif_adjust.chains[jj].residues.size();kk++) {
        if(motif_adjust.chains[jj].residues[kk].get_residue_name()==overlapped_resenames[ii]) {
          vector<double> xyz=motif_adjust.chains[jj].residues[kk].get_residue_xyz();
          Vxyz_tar.insert(Vxyz_tar.end(),xyz.begin(),xyz.end());
          check=1;
          break;
        }
      }
      if(check) break;
      //if(motif_adjust.chains[jj].residues[0].get_residue_name()==overlapped_resenames[ii]) {
        //vector<double> xyz=motif_adjust.chains[jj].residues[0].get_residue_xyz();
        //Vxyz_tar.insert(Vxyz_tar.end(),xyz.begin(),xyz.end());
        //Vxyz_interface_tar.insert(Vxyz_interface_tar.end(),xyz.begin(),xyz.end());
        //vector<double> xyz_1=motif_adjust.chains[jj].residues[1].get_residue_xyz();
        //Vxyz_interface_tar.insert(Vxyz_interface_tar.end(),xyz_1.begin(),xyz_1.end());
        //break;
      //}
      //else if(motif_adjust.chains[jj].residues.back().get_residue_name()==overlapped_resenames[ii]) {
        //size_t nres=motif_adjust.chains[jj].residues.size();
        //vector<double> xyz=motif_adjust.chains[jj].residues[nres-1].get_residue_xyz();
        //Vxyz_tar.insert(Vxyz_tar.end(),xyz.begin(),xyz.end());
        //Vxyz_interface_tar.insert(Vxyz_interface_tar.end(),xyz.begin(),xyz.end());
        //vector<double> xyz_1=motif_adjust.chains[jj].residues[nres-2].get_residue_xyz();
        //Vxyz_interface_tar.insert(Vxyz_interface_tar.end(),xyz_1.begin(),xyz_1.end());
        //break;
      //}
    }
  }

  for(size_t jj=0;jj<motif_adjust.chains.size();jj++) {
    vector<double> xyz=motif_adjust.chains[jj].get_chain_xyz();
    Vxyz_interface_tar.insert(Vxyz_interface_tar.end(),xyz.begin(),xyz.end());
  }

  vector<double> xcm_ref,xcm_adjust;
  for(size_t ii=0;ii<3;ii++) {
    xcm_ref.push_back(0.0);
    xcm_adjust.push_back(0.0);
  }
  size_t ndata=Vxyz_ref.size()/3;
  for(size_t ii=0;ii<ndata;ii++) {
    xcm_ref[0]+=Vxyz_ref[3*ii];
    xcm_ref[1]+=Vxyz_ref[3*ii+1];
    xcm_ref[2]+=Vxyz_ref[3*ii+2];
  }
  for(size_t ii=0;ii<3;ii++)
    xcm_ref[ii]/=1.0*ndata;

  ndata=Vxyz_interface_tar.size()/3;
  for(size_t ii=0;ii<ndata;ii++) {
    xcm_adjust[0]+=Vxyz_interface_tar[ii*3];
    xcm_adjust[1]+=Vxyz_interface_tar[ii*3+1];
    xcm_adjust[2]+=Vxyz_interface_tar[ii*3+2];
  }
  for(size_t ii=0;ii<3;ii++)
    xcm_adjust[ii]/=1.0*ndata;
  
  /*
  for(size_t kk=0;kk<Vmotif_reference.size();kk++)
    cout<<"Ref: "<<Vmotif_reference[kk].get_RNA_motif_info()<<endl;
  cout<<"Tar: "<<motif_target.get_RNA_motif_info()<<endl;  
  cout<<"Overlapped resenames:"<<endl;
  for(auto resname : overlapped_resenames)
    cout<<" "<<resname;
  cout<<endl;
  cout<<"Data size: "<<Vxyz_ref.size()<<" : "<<Vxyz_tar.size()<<endl;
  */

  Superpose sp_one(Vxyz_ref,Vxyz_tar);
  sp_one.calculate_rotation_matrix();
  sp_one.reset_x2(xcm_adjust);
  xcm_adjust=sp_one.rotate_x2();

  vector<double> Vxyz_adjust=motif_adjust.get_molecule_xyz();
  sp_one.reset_x2(Vxyz_adjust);
  Vxyz_adjust=sp_one.rotate_x2();
  motif_adjust.reset_molecule_xyz(Vxyz_adjust);

  vector<Atom> Vatoms_adjust=motif_adjust.get_atom_list();
  int number_of_clash_atoms=AssembleRNAMotif::calculate_number_of_clash_atoms(Vatoms_ref,Vatoms_adjust);

  double r=Measure::distance(xcm_adjust,xcm_ref);
  double delta_distance=init->assembly_delta_distance;
  double delx=(xcm_adjust[0]-xcm_ref[0])*delta_distance/r;
  double dely=(xcm_adjust[1]-xcm_ref[1])*delta_distance/r;
  double delz=(xcm_adjust[2]-xcm_ref[2])*delta_distance/r;
  if(r < MathConst::SMALL) {
    delx=(rand_one->uniform()-0.5)*2.0*delta_distance;  //[-1,1)
    dely=(rand_one->uniform()-0.5)*2.0*sqrt(delta_distance*delta_distance-delx*delx); 
    delz=sqrt(delta_distance*delta_distance-delx*delx-dely*dely);
    if(rand_one->uniform() < 0.5)
      delz*=-1.0;
  }

  int sum_move=0;
  if(flag_move) {
    while(number_of_clash_atoms>=N_CLASH_ATOMS) {
      int nmove=round(number_of_clash_atoms/N_CLASH_ATOMS)-1;
      if(nmove<1) nmove=1;
      else if(nmove>15) nmove=15;
      sum_move+=nmove;
      motif_adjust.move_molecule(delx*nmove,dely*nmove,delz*nmove);
      Vatoms_adjust.clear();
      Vatoms_adjust=motif_adjust.get_atom_list();
      number_of_clash_atoms=AssembleRNAMotif::calculate_number_of_clash_atoms(Vatoms_ref,Vatoms_adjust);
    }
  }
  
  moved_distance=sum_move*delta_distance;

  if(sum_move>=20) {
    string str="assemble_RNA_motif.cpp: Move "+to_string(moved_distance)+" Angstrom in assembling motif "+motif_target.get_RNA_motif_info();
    error->log(str);
  }
}

void AssembleRNAMotif::place_RNA_motif(const vector<RNAMotif> &Vmotif_reference, vector<RNAMotif> &Vmotif_placed) {
  double Position_Array[][3]={{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1},{1,1,0},{1,-1,0},{1,0,1},
                     {1,0,-1},{-1,1,0},{-1,-1,0},{-1,0,1},{-1,0,-1},{0,1,1},{0,1,-1},{0,-1,1},{0,-1,-1},
                     {1,1,1},{1,-1,1},{1,1,-1},{1,-1,-1},{-1,1,1},{-1,-1,1},{-1,1,-1},{-1,-1,-1}};
  size_t Num_position_array=sizeof(Position_Array)/sizeof(Position_Array[0]);
  
  vector<Atom> Vatoms_ref;
  for(size_t ii=0;ii<Vmotif_reference.size();ii++) {
    vector<Atom> Vatoms_list=Vmotif_reference[ii].get_atom_list();
    Vatoms_ref.insert(Vatoms_ref.end(),Vatoms_list.begin(),Vatoms_list.end());
  }
  vector<double> xcm_ref;
  for(size_t ii=0;ii<3;ii++)
    xcm_ref.push_back(0.0);
  for(size_t ii=0;ii<Vatoms_ref.size();ii++) {
    xcm_ref[0]+=Vatoms_ref[ii].x;
    xcm_ref[1]+=Vatoms_ref[ii].y;
    xcm_ref[2]+=Vatoms_ref[ii].z;
  }
  if(Vatoms_ref.size()>0) {
    for(size_t ii=0;ii<3;ii++)
      xcm_ref[ii]/=Vatoms_ref.size();
  }
  vector<double> Vxyz_ref;
  for(size_t ii=0;ii<Vmotif_reference.size();ii++) {
    int flag_break=0;
    for(size_t jj=0;jj<Vmotif_reference[ii].chains.size();jj++) {
      if(Vmotif_reference[ii].chains[jj].residues.size()>=2) {
        Residue res_1=Vmotif_reference[ii].chains[jj].residues[0];
        Residue res_2=Vmotif_reference[ii].chains[jj].residues[1];
        res_1.simplify_RNA_base();
        res_2.simplify_RNA_base();
        vector<double> xyz1=res_1.get_RNA_residue_xyz();
        vector<double> xyz2=res_2.get_RNA_residue_xyz();
        Vxyz_ref.insert(Vxyz_ref.end(),xyz1.begin(),xyz1.end());
        Vxyz_ref.insert(Vxyz_ref.end(),xyz2.begin(),xyz2.end());
        flag_break=1;
        break;
      }
    }
    if(flag_break) break;
  }
  if(Vxyz_ref.size()==0) {
    string str="assemble_RNA_motif.cpp: Chains are too short in reference motifs!";
    error->one(str);
  }

  for(size_t ii=0;ii<Vmotif_placed.size();ii++) {
    motif_adjust=Vmotif_placed[ii];
    vector<double> Vxyz_tar;
    for(size_t jj=0;jj<motif_adjust.chains.size();jj++) {
      if(motif_adjust.chains[jj].residues.size()>=2) {
        Residue res_1=motif_adjust.chains[jj].residues[0];
        Residue res_2=motif_adjust.chains[jj].residues[1];
        res_1.simplify_RNA_base();
        res_2.simplify_RNA_base();
        vector<double> xyz1=res_1.get_RNA_residue_xyz();
        vector<double> xyz2=res_2.get_RNA_residue_xyz();
        Vxyz_tar.insert(Vxyz_tar.end(),xyz1.begin(),xyz1.end());
        Vxyz_tar.insert(Vxyz_tar.end(),xyz2.begin(),xyz2.end());
        break;
      }
    }
    if(Vxyz_tar.size()==0) {
      string str="assemble_RNA_motif.cpp: Chains are too short in motif "+motif_adjust.get_RNA_motif_info();
      error->one(str);
    }

    Superpose sp_one(Vxyz_ref,Vxyz_tar);
    sp_one.calculate_rotation_matrix();
    vector<double> Vxyz_adjust=motif_adjust.get_molecule_xyz();
    sp_one.reset_x2(Vxyz_adjust);
    Vxyz_adjust=sp_one.rotate_x2();
    motif_adjust.reset_molecule_xyz(Vxyz_adjust);

    vector<double> xcm_adjust;
    for(size_t kk=0;kk<3;kk++)
      xcm_adjust.push_back(0.0);
    size_t natom_adjust=Vxyz_adjust.size()/3;
    for(size_t jj=0;jj<natom_adjust;jj++) {
      xcm_adjust[0]+=Vxyz_adjust[3*jj];
      xcm_adjust[1]+=Vxyz_adjust[3*jj+1];
      xcm_adjust[2]+=Vxyz_adjust[3*jj+2];
    }
    for(size_t kk=0;kk<3;kk++)
      xcm_adjust[kk]/=natom_adjust;
    vector<double> xyz_move;
    for(size_t kk=0;kk<3;kk++) {
      xyz_move.push_back(xcm_ref[kk]-xcm_adjust[kk]);
    }
    motif_adjust.move_molecule(xyz_move[0],xyz_move[1],xyz_move[2]);

    int index_i=ii%Num_position_array;
    xyz_move.clear();
    for(size_t kk=0;kk<3;kk++) {
      xyz_move.push_back(Position_Array[index_i][kk]);
    }
    double r=MathExtra::length3(xyz_move);
    for(size_t kk=0;kk<3;kk++) {
      xyz_move[kk]=xyz_move[kk]*init->assembly_delta_distance/r;
    }
    
    int number_of_clash_atoms=2*N_CLASH_ATOMS;
    int sum_move=0;
    while(number_of_clash_atoms>=N_CLASH_ATOMS) {
      int nmove=round(number_of_clash_atoms/N_CLASH_ATOMS)-1;
      if(nmove<1) nmove=1;
      else if(nmove>15) nmove=15;
      sum_move+=nmove;
      motif_adjust.move_molecule(xyz_move[0]*nmove,xyz_move[1]*nmove,xyz_move[2]*nmove);
      vector<Atom> Vatoms_adjust=motif_adjust.get_atom_list();
      number_of_clash_atoms=AssembleRNAMotif::calculate_number_of_clash_atoms(Vatoms_ref,Vatoms_adjust);
    }
  
    moved_distance=sum_move*init->assembly_delta_distance;

    //if(sum_move>0) {
      //string str="assemble_RNA_motif.cpp: Move "+to_string(moved_distance)+" Angstrom in placing motif "+motif_adjust.get_RNA_motif_info();
      //error->log(str);
    //}

    Vmotif_placed[ii]=motif_adjust;

    vector<Atom> Vatoms_list=motif_adjust.get_atom_list();
    Vatoms_ref.insert(Vatoms_ref.end(),Vatoms_list.begin(),Vatoms_list.end());
  }
}

void AssembleRNAMotif::place_RNA_motif(const vector<RNAMotif> &Vmotif_reference, vector<vector<RNAMotif> > &VVmotif_placed) {
  double Position_Array[][3]={{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1},{1,1,0},{1,-1,0},{1,0,1},
                     {1,0,-1},{-1,1,0},{-1,-1,0},{-1,0,1},{-1,0,-1},{0,1,1},{0,1,-1},{0,-1,1},{0,-1,-1},
                     {1,1,1},{1,-1,1},{1,1,-1},{1,-1,-1},{-1,1,1},{-1,-1,1},{-1,1,-1},{-1,-1,-1}};
  size_t Num_position_array=sizeof(Position_Array)/sizeof(Position_Array[0]);
  
  vector<Atom> Vatoms_ref;
  for(size_t ii=0;ii<Vmotif_reference.size();ii++) {
    vector<Atom> Vatoms_list=Vmotif_reference[ii].get_atom_list();
    Vatoms_ref.insert(Vatoms_ref.end(),Vatoms_list.begin(),Vatoms_list.end());
  }
  vector<double> xcm_ref;
  for(size_t ii=0;ii<3;ii++)
    xcm_ref.push_back(0.0);
  for(size_t ii=0;ii<Vatoms_ref.size();ii++) {
    xcm_ref[0]+=Vatoms_ref[ii].x;
    xcm_ref[1]+=Vatoms_ref[ii].y;
    xcm_ref[2]+=Vatoms_ref[ii].z;
  }
  if(Vatoms_ref.size()>0) {
    for(size_t ii=0;ii<3;ii++)
      xcm_ref[ii]/=Vatoms_ref.size();
  }

  for(size_t ii=0;ii<VVmotif_placed.size();ii++) {
    vector<Atom> Vatoms_adjust;
    for(size_t jj=0;jj<VVmotif_placed[ii].size();jj++) {
      vector<Atom> Vatoms_list=VVmotif_placed[ii][jj].get_atom_list();
      Vatoms_adjust.insert(Vatoms_adjust.end(),Vatoms_list.begin(),Vatoms_list.end());
    }
    vector<double> xcm_adjust;
    for(size_t kk=0;kk<3;kk++)
      xcm_adjust.push_back(0.0);
    size_t natom_adjust=Vatoms_adjust.size();
    for(size_t jj=0;jj<natom_adjust;jj++) {
      xcm_adjust[0]+=Vatoms_adjust[jj].x;
      xcm_adjust[1]+=Vatoms_adjust[jj].y;
      xcm_adjust[2]+=Vatoms_adjust[jj].z;
    }
    if(natom_adjust>0) {
      for(size_t kk=0;kk<3;kk++)
        xcm_adjust[kk]/=natom_adjust;
    }
    vector<double> xyz_move_sum;
    for(size_t kk=0;kk<3;kk++) {
      xyz_move_sum.push_back(xcm_ref[kk]-xcm_adjust[kk]);
    }
    for(size_t jj=0;jj<natom_adjust;jj++) {
      Vatoms_adjust[jj].move_atom(xyz_move_sum[0],xyz_move_sum[1],xyz_move_sum[2]);
    }

    int index_i=ii%Num_position_array;
    vector<double> xyz_move;
    for(size_t kk=0;kk<3;kk++) {
      xyz_move.push_back(Position_Array[index_i][kk]);
    }
    double r=MathExtra::length3(xyz_move);
    for(size_t kk=0;kk<3;kk++) {
      xyz_move[kk]=xyz_move[kk]*init->assembly_delta_distance/r;
    }
    
    int number_of_clash_atoms=2*N_CLASH_ATOMS;
    int sum_move=0;
    while(number_of_clash_atoms>=N_CLASH_ATOMS) {
      int nmove=round(number_of_clash_atoms/N_CLASH_ATOMS)-1;
      if(nmove<1) nmove=1;
      else if(nmove>15) nmove=15;
      sum_move+=nmove;
      for(size_t jj=0;jj<natom_adjust;jj++) {
        Vatoms_adjust[jj].move_atom(xyz_move[0]*nmove,xyz_move[1]*nmove,xyz_move[2]*nmove);
      }
      number_of_clash_atoms=AssembleRNAMotif::calculate_number_of_clash_atoms(Vatoms_ref,Vatoms_adjust);
    }
  
    for(int kk=0;kk<3;kk++)
      xyz_move_sum[kk]+=sum_move*xyz_move[kk];

    for(size_t jj=0;jj<VVmotif_placed[ii].size();jj++) {
      VVmotif_placed[ii][jj].move_molecule(xyz_move_sum[0],xyz_move_sum[1],xyz_move_sum[2]);
    }

    Vatoms_ref.insert(Vatoms_ref.end(),Vatoms_adjust.begin(),Vatoms_adjust.end());
  }
}

int AssembleRNAMotif::calculate_number_of_clash_atoms(const vector<Atom> &Vatoms_ref, const vector<Atom> &Vatoms_adjust) {
  int number_of_clash_atoms=0;
  for(size_t ii=0;ii<Vatoms_ref.size();ii++) {
    string atomname_i=Vatoms_ref[ii].get_atom_name();
    for(size_t jj=0;jj<Vatoms_adjust.size();jj++) {
      string atomname_j=Vatoms_adjust[jj].get_atom_name();
      if(atomname_i==atomname_j) continue;
      double distance=Measure::distance(Vatoms_ref[ii],Vatoms_adjust[jj]);
      if(distance<ATOM_DISTANCE_MIN)
        number_of_clash_atoms++;
    }
  }

  return number_of_clash_atoms;
}

void AssembleRNAMotif::clear() {
  motif_adjust.clear();
  overlapped_resenames.clear();
  overlapped_resenames.shrink_to_fit();
  moved_distance=0.0;
}