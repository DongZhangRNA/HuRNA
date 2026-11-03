#include "basic.h"
#include "readfile.h"
#include "initialization.h"
#include "error.h"
#include "superpose.h"
#include "RNA_molecule.h"
#include "RNA_basepair.h"

using namespace CBS_NS;

RNAMolecule::RNAMolecule(Cobmms *cbs) : Molecule(cbs) {}

RNAMolecule::RNAMolecule(Cobmms *cbs, string seq_one) : Molecule(cbs) {
  RNAMolecule::set_RNA_molecule(seq_one);
}

RNAMolecule::RNAMolecule(const RNAMolecule &other) : Molecule(other) {
  sequences=other.sequences;
  basepairs_WC=other.basepairs_WC;
  basepairs_nWC=other.basepairs_nWC;
  basestacks=other.basestacks;
}

RNAMolecule& RNAMolecule::operator=(const RNAMolecule &other) {
  if(this != &other) {
    Molecule::operator=(other);
    sequences=other.sequences;
    basepairs_WC=other.basepairs_WC;
    basepairs_nWC=other.basepairs_nWC;
    basestacks=other.basestacks;
  }

  return *this;
}

RNAMolecule::~RNAMolecule() {
  sequences.clear();
  sequences.shrink_to_fit();
  chain_resnames.clear();
  chain_resnames.shrink_to_fit();
  chains.clear();
  chains.shrink_to_fit();
  chain_ids.clear();
  chain_ids.shrink_to_fit();
  chain_labels.clear();
  chain_labels.shrink_to_fit();
  Vresnames.clear();
  Vresnames.shrink_to_fit();
  Vchainid.clear();
  Vchainid.shrink_to_fit();
  basepairs_WC.clear();
  basepairs_WC.shrink_to_fit();
  basepairs_nWC.clear();
  basepairs_nWC.shrink_to_fit();
  basestacks.clear();
  basestacks.shrink_to_fit();
}

void RNAMolecule::set_RNA_molecule(const RNAMolecule &rna_one) {
  sequences=rna_one.sequences;
  chain_resnames=rna_one.chain_resnames;
  chains=rna_one.chains;
  chain_ids=rna_one.chain_ids;
  chain_labels=rna_one.chain_labels;
  Vresnames=rna_one.Vresnames;
  Vchainid=rna_one.Vchainid;
  basepairs_WC=rna_one.basepairs_WC;
  basepairs_nWC=rna_one.basepairs_nWC;
  basestacks=rna_one.basestacks;
}

void RNAMolecule::set_RNA_molecule(string seq_one, int chainid_offset) {
  sequences.clear();
  sequences=RNAMolecule::read_sequence(seq_one);

  vector<vector<string> > Vresname_one;
  int nres=0;
  for(int ii=0;ii<(int)sequences.size();ii++) {
    int cid=(ii+chainid_offset)%((int)init->CHAIN_ID_LIST.size());
    int nround=(ii+chainid_offset)/((int)init->CHAIN_ID_LIST.size());
    string chainid_one=init->CHAIN_ID_LIST[cid];
    if(nround>0) chainid_one+=init->CHAIN_ID_LIST[nround-1];
    vector<string> Vtmp=RNAMolecule::get_resnames_from_sequence(sequences[ii]);
    for(size_t jj=0;jj<Vtmp.size();jj++) {
      nres++;
      string resname_one=chainid_one+"/"+to_string(nres)+"/"+Vtmp[jj];
      Vtmp[jj]=resname_one;
    }
    Vresname_one.push_back(Vtmp);
  }
  
  Molecule::set_molecule(Vresname_one,Chain::rna);
}

void RNAMolecule::set_RNA_molecule(const vector<Atom> &Vatoms) {
  sequences.clear();
  Molecule::clear();

  vector<Residue> Vnts;
  vector<Atom> nt_atoms;
  string nt_name0=Vatoms[0].chainid+"/"+to_string(Vatoms[0].resid);
  if(!isspace(Vatoms[0].iCode))
    nt_name0+=to_string((char)Vatoms[0].iCode);
  nt_name0+="/"+Vatoms[0].resname;
  for(size_t ii=0;ii<Vatoms.size();ii++) {
    string nt_name=Vatoms[ii].chainid+"/"+to_string(Vatoms[ii].resid);
    if(!isspace(Vatoms[ii].iCode))
      nt_name+=to_string((char)Vatoms[ii].iCode);
    nt_name+="/"+Vatoms[ii].resname;
    if(nt_name!=nt_name0) {
      Residue nt(cbs);
      nt.set_residue(nt_atoms);
      if(nt.check_atom("C4'") && nt.check_atom("C3'") && (nt.check_atom("O3'") || nt.check_atom("C5'"))) {
        nt.recovery_RNA_backbone_atoms();
        nt.recovery_RNA_base();
        Vnts.push_back(nt);
      }
      nt_name0=nt_name;
      nt_atoms.clear();
    }
    nt_atoms.push_back(Vatoms[ii]);
  }
  if(nt_atoms.size()>4) {
    Residue nt(cbs);
    nt.set_residue(nt_atoms);
    if(nt.check_atom("C4'") && nt.check_atom("C3'") && (nt.check_atom("O3'") || nt.check_atom("C5'"))) {
      nt.recovery_RNA_backbone_atoms();
      nt.recovery_RNA_base();
      Vnts.push_back(nt);
    }
  }

  if(Vnts.size()==0) return;

  vector<Residue> chain_nts;
  chain_nts.push_back(Vnts[0]);
  for(size_t ii=1;ii<Vnts.size();ii++) {
    int chain_break=0;
    if(!(Vnts[ii-1].check_atom("O3'")) && !(Vnts[ii].check_atom("P"))) chain_break=1;
    else {
      Atom atom_P=Vnts[ii]["P"];
      Atom atom_O3p=Vnts[ii-1]["O3'"];
      double dis_OP=atom_P.distance(atom_O3p);
      if(dis_OP < init->RNA_distance_O3p_P_min || dis_OP > init->RNA_distance_O3p_P_max) chain_break=1;
    }
    if(Vnts[ii].chainid!=Vnts[ii-1].chainid || (init->RNA_check_chain_break && chain_break)) {
      if(chain_nts.size()>0) {
        Chain chain(cbs);
        chain.set_chain(chain_nts,Chain::rna);
        Molecule::insert_chain(chain);
      }
      chain_nts.clear();
    }
    chain_nts.push_back(Vnts[ii]);
  }
  if(chain_nts.size()>0) {
    Chain chain(cbs);
    chain.set_chain(chain_nts,Chain::rna);
    Molecule::insert_chain(chain);
  }

  for(size_t ii=0;ii<chains.size();ii++) {
    string seq=chains[ii].get_sequence();
    sequences.push_back(seq);
  }
}

string RNAMolecule::get_sequences() const{
  string seq=sequences[0];
  for(size_t ii=1;ii<sequences.size();ii++) {
    seq+="&"+sequences[ii];
  }

  return seq;
}

void RNAMolecule::replace_sequence(string &seq_one) {
  vector<string> Vseq=RNAMolecule::read_sequence(seq_one);

  size_t ni=sequences.size();
  if(ni!=Vseq.size()) {
    string str="RNA_molecule.cpp: the new sequences does not match the number of chains in the replaced molecule --> ";
    str+=to_string((int) ni)+"/"+to_string((int) Vseq.size());
    error->one(str);
  }

  for(size_t ii=0;ii<ni;ii++) {
    if(Vseq[ii]=="N") {
      continue;
    }
    else {
      size_t nres_old=RNAMolecule::get_resnames_from_sequence(sequences[ii]).size();
      if(Vseq[ii].size()!=nres_old) {
        string str="RNA_molecule.cpp: the length of new sequence does not match the length of the replaced molecule --> ";
        str+="New seq: "+Vseq[ii]+"/n"+"Old seq:"+sequences[ii];
        error->one(str);
      }
      sequences[ii]=Vseq[ii];
    }
  }

  vector<vector<string> > Vseq_one;
  for(size_t ii=0;ii<sequences.size();ii++) {
    vector<string> Vtmp=RNAMolecule::get_resnames_from_sequence(sequences[ii]);
    Vseq_one.push_back(Vtmp);
  }

  Molecule::replace_sequence(Vseq_one);
}

string RNAMolecule::extract_secondary_structure() const{
  Molecule mol_one=Molecule::copy_molecule();
  mol_one.standardize_molecule();
  
  string fprefix=init->WORK_DIR+"dssr";
  string fstr=fprefix+"-input.pdb";
  string fout=fprefix+"-output.dat";
  ofstream outF (fstr.c_str());
  outF<<mol_one;
  outF.close(); outF.clear();

  string cmd=init->DSSR+" -i="+fstr+" -o="+fout+" --prefix="+fprefix;
  cbs->execute_system_command(cmd);

  string fdbn=fprefix+"-2ndstrs.dbn";
  if(!is_file_exist(fdbn)) {
    string str="RNA_molecule.cpp: Failed to extract secondary structure!";
    error->warning(str);

    string secstr="";
    for(size_t ii=0;ii<Molecule::number_of_residues();ii++)
      secstr+=".";
    
    return secstr;
  }
  
  vector<string> Vstr=read_file(fdbn);
  cmd="rm "+fprefix+"-*";
  cbs->execute_system_command(cmd);

  if(Vstr[1].size()!=RNAMolecule::get_sequences().size()) {
    string str="RNA_molecule.cpp: Sequence given by DSSR is not match --> \n";
    str+=RNAMolecule::get_sequences()+"\n"+Vstr[1];
    error->warning(str);
  }

  return Vstr[2];
}

void RNAMolecule::extract_interactions() {
  basepairs_WC.clear();
  basepairs_nWC.clear();
  basestacks.clear();

  string fprefix=init->WORK_DIR+"dssr";
  string fstr=fprefix+"-input.pdb";
  ofstream outF (fstr.c_str());
  int index=0;
  for(size_t ii=0;ii<chains.size();ii++) {
    index=chains[ii].set_atomid(index)+1;
    outF<<chains[ii];
  }
  outF.close(); outF.clear();

  string fout=fprefix+"-output.dat";
  string cmd=init->DSSR+" -i="+fstr+" -o="+fout+" --prefix="+fprefix;
  cbs->execute_system_command(cmd);

  string fdbn=fprefix+"-2ndstrs.dbn";
  if(!is_file_exist(fdbn)) {
    string str="RNA_molecule.cpp: Failed to extract secondary structure!";
    error->warning(str);
    return;
  }

  vector<string> Vstr=read_file(fout);
  vector<string> Vsstr=read_file(fdbn);
  cmd="rm "+fprefix+"-*";
  cbs->execute_system_command(cmd);

  size_t ni=0;
  while(ni<Vstr.size()) {
    if(Vstr[ni].find("List of")==string::npos) {
      ni++;
      continue;
    }
    vector<string> keys=split_string(Vstr[ni]);
    if(keys[0]!="List" && keys[1]!="of") {
      ni++;
      continue;
    }
    if(keys[3]=="base" && (keys[4]=="pairs"||keys[4]=="pair")) {
      int nbp=stoi(keys[2]);
      ni++;
      for(int jj=1;jj<=nbp;jj++) {
        ni++;
        vector<string> stoken=split_string(Vstr[ni]);
        if(stoken.size()!=8) {
          string str="RNA_molecule.cpp: Wrong basepair line info --> "+Vstr[ni];
          error->one(str);
        }
        if(jj!=stoi(stoken[0])) {
          string str="RNA_molecule.cpp: Ids for basepair do not match --> "+to_string(jj)+"/"+stoken[0];
          error->message(str);
        }
        string nt1=stoken[1];
        string nt2=stoken[2];
        string bp=stoken[3];
        string LW=stoken[6];
        string res_i=RNAMolecule::change_dssr_residue_name(nt1);
        int index_i=Molecule::get_residue_index(res_i);
        string res_j=RNAMolecule::change_dssr_residue_name(nt2);
        int index_j=Molecule::get_residue_index(res_j);
        string type=to_string((char)bp[0])+to_string((char)bp[2])+"_"+LW;
        InteractionPair bp_i;
        bp_i.residue_i=res_i;
        bp_i.residue_j=res_j;
        bp_i.type=type;
        if(index_i>index_j) RNAMolecule::swap_interactionpair_type(bp_i);
        if(RNAMolecule::is_WC_basepair(type)) basepairs_WC.push_back(bp_i);
        else basepairs_nWC.push_back(bp_i);
      }
    }
    else if(keys[3]=="stems" || keys[3]=="stem") {
      int nstem=stoi(keys[2]);
      ni+=5;
      for(int jj=0;jj<nstem;jj++) {
        ni++;
        size_t found=Vstr[ni].find("bps=");
        if(found==string::npos) {
          string str="RNA_molecule.cpp: Wrong stem head line --> "+Vstr[ni];
          error->one(str);
        }
        string sline=Vstr[ni].substr(found);
        vector<string> stoken=split_string(sline,"=");
        int nbps=stoi(stoken[1]);
        ni+=5;
        vector<string> Vnt1,Vnt2;
        for(int kk=1;kk<=nbps;kk++) {
          stoken=split_string(Vstr[ni]);
          if(stoken.size()!=8) {
            string str="RNA_molecule.cpp: Wrong stem line info --> "+Vstr[ni];
            error->one(str);
          }
          if(kk!=stoi(stoken[0])) {
            string str="RNA_molecule.cpp: Ids for stem basepair do not match --> "+to_string(kk)+"/"+stoken[0];
            error->message(str);
          }
          Vnt1.push_back(stoken[1]);
          Vnt2.push_back(stoken[2]);
          ni++;
        }
        for(size_t kk=0;kk<Vnt1.size()-1;kk++) {
          string res_i=RNAMolecule::change_dssr_residue_name(Vnt1[kk]);
          int index_i=Molecule::get_residue_index(res_i);
          string res_j=RNAMolecule::change_dssr_residue_name(Vnt1[kk+1]);
          int index_j=Molecule::get_residue_index(res_j);
          stoken=split_string(res_i,"/");
          string type=stoken[2];
          stoken=split_string(res_j,"/");
          type+=stoken[2];
          InteractionPair bs_i;
          bs_i.residue_i=res_i;
          bs_i.residue_j=res_j;
          bs_i.type=type;
          if(index_i>index_j) RNAMolecule::swap_interactionpair_type(bs_i);
          basestacks.push_back(bs_i);
        }
        for(int kk=(int)Vnt2.size()-1;kk>0;kk--) {
          string res_i=RNAMolecule::change_dssr_residue_name(Vnt2[kk]);
          int index_i=Molecule::get_residue_index(res_i);
          string res_j=RNAMolecule::change_dssr_residue_name(Vnt2[kk-1]);
          int index_j=Molecule::get_residue_index(res_j);
          stoken=split_string(res_i,"/");
          string type=stoken[2];
          stoken=split_string(res_j,"/");
          type+=stoken[2];
          InteractionPair bs_i;
          bs_i.residue_i=res_i;
          bs_i.residue_j=res_j;
          bs_i.type=type;
          if(index_i>index_j) RNAMolecule::swap_interactionpair_type(bs_i);
          basestacks.push_back(bs_i);
        }
      }
    }
    else if(keys[3]=="stacks" || keys[3]=="stack") {
      int nstack=stoi(keys[2]);
      ni+=4;
      for(int jj=1;jj<=nstack;jj++) {
        vector<string> stoken=split_string(Vstr[ni],",=");
        if(jj!=stoi(stoken[0])) {
          string str="RNA_molecule.cpp: Ids for stack do not match --> "+to_string(jj)+"/"+stoken[0];
          error->message(str);
        }
        int nnt=stoi(stoken[2]);
        if((int)stoken.size()!=nnt+4) {
          string str="RNA_molecule.cpp: Wrong stack line info --> "+Vstr[ni];
          error->one(str);
        }
        vector<string> Vnt;
        for(size_t kk=4;kk<stoken.size();kk++)
          Vnt.push_back(stoken[kk]);
        ni++;
        for(size_t kk=0;kk<Vnt.size()-1;kk++) {
          string res_i=RNAMolecule::change_dssr_residue_name(Vnt[kk]);
          int index_i=Molecule::get_residue_index(res_i);
          string res_j=RNAMolecule::change_dssr_residue_name(Vnt[kk+1]);
          int index_j=Molecule::get_residue_index(res_j);
          stoken=split_string(res_i,"/");
          string type=stoken[2];
          stoken=split_string(res_j,"/");
          type+=stoken[2];
          InteractionPair bs_i;
          bs_i.residue_i=res_i;
          bs_i.residue_j=res_j;
          bs_i.type=type;
          if(index_i>index_j) RNAMolecule::swap_interactionpair_type(bs_i);
          basestacks.push_back(bs_i);
        }
      }
    }
    else
      ni++;
  }
}

string RNAMolecule::extract_secondary_structure_MCAnnotate() const{
  Molecule mol_one=Molecule::copy_molecule();
  mol_one.standardize_molecule();
  
  string fprefix=init->WORK_DIR+"mcannotate";
  string fstr=fprefix+"-input.pdb";
  ofstream outF (fstr.c_str());
  int index=0;
  for(size_t ii=0;ii<mol_one.chains.size();ii++) {
    index=mol_one.chains[ii].set_atomid(index)+1;
    outF<<mol_one.chains[ii];
  }
  outF.close(); outF.clear();

  string fout=fprefix+"-output.dat";
  string cmd=init->MCAnnotate+" "+fstr+" > "+fout;
  cbs->execute_system_command(cmd);

  vector<string> Vlines=read_file(fout);
  cmd="rm "+fprefix+"-*";
  cbs->execute_system_command(cmd);

  enum State { STATE_OUT = 0, STATE_RESIDUE = 1, STATE_PAIR = 2, STATE_STACK = 3 };
  static const std::regex pattern_residue(
    R"(^([\w']+) : (\w+))",
    std::regex_constants::optimize
  );
  static const std::regex pattern_pair(
    R"(^([\w']+)-([\w']+) : (\w+)-(\w+) ([\w']+)/([\w']+).*pairing( (parallel|antiparallel) (cis|trans)))",
    std::regex_constants::optimize
  );
  
  int model_count = 0;
  int state = STATE_OUT;
  vector<string> Vres_mca,Vseq;
  map<string,string> Mbp_mca;
  for(size_t ii=0;ii<Vlines.size();ii++) {
    string sline=trim(Vlines[ii]);
    // Only first model used
    if(line_starts_with(sline, "Residue conformations")) {
      if(model_count == 0) {
        state = STATE_RESIDUE;
        ++model_count;
        continue;
      } 
      else {
        break;
      }
    }

    if(line_starts_with(sline, "Base-pairs")) {
      state = STATE_PAIR;
      continue;
    }

    if(line_starts_with(sline, "Adjacent stackings") || line_starts_with(sline, "Non-Adjacent stackings")) {
      state = STATE_STACK;
      continue;
    }

    if(line_ends_with(sline, "----------")) {
      state = STATE_OUT;
      continue;
    }

    if(state == STATE_RESIDUE) {
      std::smatch sm;
      if(std::regex_search(sline, sm, pattern_residue)) {
        Vres_mca.push_back(sm[1]);
        Vseq.push_back(sm[2]);
      }
    }

    if(state == STATE_PAIR) {
      std::smatch sm;
      if(std::regex_search(sline, sm, pattern_pair)) {
        char side_i=(char)std::toupper((unsigned char)sm[5].str()[0]);
        char side_j=(char)std::toupper((unsigned char)sm[6].str()[0]);
        if (!((side_i == 'W' || side_i == 'H' || side_i == 'S') &&
              (side_j == 'W' || side_j == 'H' || side_j == 'S'))) {
          continue;
        }
        string nt_i=sm[1].str();
        string nt_j=sm[2].str();
        string type=sm[3].str()+sm[4].str();
        string bp_orientation=sm[9].str();
        for(size_t kk=0;kk<bp_orientation.size();kk++) 
          bp_orientation[kk]=std::tolower(bp_orientation[kk]);
        if(bp_orientation=="cis")
          type+="_c";
        else if(bp_orientation=="trans")
          type+="_t";
        type+=to_string((char)side_i)+to_string((char)side_j);
        if(RNAMolecule::is_WC_basepair(type)) 
          Mbp_mca.emplace(nt_i,nt_j);
      }
    }
  }

  map<string, string> Mres_mca_std;
  if(Vres_mca.size()!=Molecule::Vresnames.size()) {
    string str="RNA_molecule.cpp: Number of residues identified by MC-Annotate is not matched with the original number! --> ";
    str+=to_string(Vres_mca.size())+" : "+to_string(Molecule::Vresnames.size());
    error->one(str);
  }
  for(size_t ii=0;ii<Vres_mca.size();ii++) {
    string res_mac=Vres_mca[ii];
    string seq_mca=Vseq[ii];
    string res_std=Molecule::Vresnames[ii];
    vector<string> stoken=split_string(res_std,"/");
    if(stoken[2]!=seq_mca) {
      string str="RNA_molecule.cpp: Residue name identified by MC-Annotate is not matched with the original one! --> ";
      str+=seq_mca+"("+res_mac+") : "+stoken[2]+"("+res_std+")";
      error->one(str);
    }
    Mres_mca_std.emplace(res_mac,res_std);
  }
  Vres_mca.clear();
  Vseq.clear();

  map<string,string> Mbps;
  for(auto& it : Mbp_mca) {
    string nt_i,nt_j;
    if(Mres_mca_std.count(it.first))
      nt_i=Mres_mca_std[it.first];
    else {
      string str="RNA_molecule.cpp: Cannot find standard residue name for MC-Annotate resiude! --> "+it.first;
      error->one(str);
    }
    if(Mres_mca_std.count(it.second))
      nt_j=Mres_mca_std[it.second];
    else {
      string str="RNA_molecule.cpp: Cannot find standard residue name for MC-Annotate resiude! --> "+it.second;
      error->one(str);
    }
    Mbps.emplace(nt_i,nt_j);
  }

  if(Mbps.size()==0) {
    string str="RNA_molecule.cpp: No base pairs found!";
    error->warning(str);

    string secstr="";
    for(size_t ii=0;ii<Molecule::number_of_residues();ii++)
      secstr+=".";
    
    return secstr;
  }
  
  RNABasePair rna_bp_one(cbs);
  vector<vector<string> > Vres_one=mol_one.get_residue_name();
  rna_bp_one.set_RNA_basepairs(Vres_one,Mbps);

  return rna_bp_one.secstr_out;
}

void RNAMolecule::extract_interactions_MCAnnotate() {
  basepairs_WC.clear();
  basepairs_nWC.clear();
  basestacks.clear();

  string fprefix=init->WORK_DIR+"mcannotate";
  string fstr=fprefix+"-input.pdb";
  ofstream outF (fstr.c_str());
  int index=0;
  for(size_t ii=0;ii<chains.size();ii++) {
    index=chains[ii].set_atomid(index)+1;
    outF<<chains[ii];
  }
  outF.close(); outF.clear();

  string fout=fprefix+"-output.dat";
  string cmd=init->MCAnnotate+" "+fstr+" > "+fout;
  cbs->execute_system_command(cmd);

  vector<string> Vlines=read_file(fout);
  cmd="rm "+fprefix+"-*";
  cbs->execute_system_command(cmd);

  enum State { STATE_OUT = 0, STATE_RESIDUE = 1, STATE_PAIR = 2, STATE_STACK = 3 };
  static const std::regex pattern_residue(
    R"(^([\w']+) : (\w+))",
    std::regex_constants::optimize
  );
  static const std::regex pattern_pair(
    R"(^([\w']+)-([\w']+) : (\w+)-(\w+) ([\w']+)/([\w']+).*pairing( (parallel|antiparallel) (cis|trans)))",
    std::regex_constants::optimize
  );
  static const std::regex pattern_stack(
    R"(^([\w']+)-([\w']+) :.*(inward|upward|downward|outward))",
    std::regex_constants::optimize
  );
  
  int model_count = 0;
  int state = STATE_OUT;
  vector<string> Vres,Vseq;
  for(size_t ii=0;ii<Vlines.size();ii++) {
    string sline=trim(Vlines[ii]);
    // Only first model used
    if(line_starts_with(sline, "Residue conformations")) {
      if(model_count == 0) {
        state = STATE_RESIDUE;
        ++model_count;
        continue;
      } 
      else {
        break;
      }
    }

    if(line_starts_with(sline, "Base-pairs")) {
      state = STATE_PAIR;
      continue;
    }

    if(line_starts_with(sline, "Adjacent stackings") || line_starts_with(sline, "Non-Adjacent stackings")) {
      state = STATE_STACK;
      continue;
    }

    if(line_ends_with(sline, "----------")) {
      state = STATE_OUT;
      continue;
    }

    if(state == STATE_RESIDUE) {
      std::smatch sm;
      if(std::regex_search(sline, sm, pattern_residue)) {
        Vres.push_back(sm[1]);
        Vseq.push_back(sm[2]);
      }
    }

    if(state == STATE_PAIR) {
      std::smatch sm;
      if(std::regex_search(sline, sm, pattern_pair)) {
        char side_i=(char)std::toupper((unsigned char)sm[5].str()[0]);
        char side_j=(char)std::toupper((unsigned char)sm[6].str()[0]);
        if (!((side_i == 'W' || side_i == 'H' || side_i == 'S') &&
              (side_j == 'W' || side_j == 'H' || side_j == 'S'))) {
          continue;
        }
        InteractionPair bp_i;
        bp_i.residue_i=sm[1].str();
        bp_i.residue_j=sm[2].str();
        string type=sm[3].str()+sm[4].str();
        string bp_orientation=sm[9].str();
        for(size_t kk=0;kk<bp_orientation.size();kk++) 
          bp_orientation[kk]=std::tolower(bp_orientation[kk]);
        if(bp_orientation=="cis")
          type+="_c";
        else if(bp_orientation=="trans")
          type+="_t";
        type+=to_string((char)side_i)+to_string((char)side_j);
        bp_i.type=type;
        if(RNAMolecule::is_WC_basepair(type)) basepairs_WC.push_back(bp_i);
        else basepairs_nWC.push_back(bp_i);
      }
    }

    if(state == STATE_STACK) {
      std::smatch sm;
      if(std::regex_search(sline, sm, pattern_stack)) {
        InteractionPair bs_i;
        bs_i.residue_i=sm[1].str();
        bs_i.residue_j=sm[2].str();
        bs_i.type=sm[3].str();
        basestacks.push_back(bs_i);
      }
    }
  }

  map<string, string> Mres_mca_std;
  if(Vres.size()!=Molecule::Vresnames.size()) {
    string str="RNA_molecule.cpp: Number of residues identified by MC-Annotate is not matched with the original number! --> ";
    str+=to_string(Vres.size())+" : "+to_string(Molecule::Vresnames.size());
    error->one(str);
  }
  for(size_t ii=0;ii<Vres.size();ii++) {
    string res_mac=Vres[ii];
    string seq_mca=Vseq[ii];
    string res_std=Molecule::Vresnames[ii];
    vector<string> stoken=split_string(res_std,"/");
    if(stoken[2]!=seq_mca) {
      string str="RNA_molecule.cpp: Residue name identified by MC-Annotate is not matched with the original one! --> ";
      str+=seq_mca+"("+res_mac+") : "+stoken[2]+"("+res_std+")";
      error->one(str);
    }
    Mres_mca_std.emplace(res_mac,res_std);
  }
  Vres.clear();
  Vseq.clear();

  for(size_t ii=0;ii<basepairs_WC.size();ii++) {
    InteractionPair IP_one;
    string res_mac_i=basepairs_WC[ii].residue_i;
    string res_mac_j=basepairs_WC[ii].residue_j;
    if(Mres_mca_std.count(res_mac_i))
      IP_one.residue_i=Mres_mca_std[res_mac_i];
    else {
      string str="RNA_molecule.cpp: Cannot find standard residue name for MC-Annotate resiude! --> "+res_mac_i;
      error->one(str);
    }
    if(Mres_mca_std.count(res_mac_j))
      IP_one.residue_j=Mres_mca_std[res_mac_j];
    else {
      string str="RNA_molecule.cpp: Cannot find standard residue name for MC-Annotate resiude! --> "+res_mac_j;
      error->one(str);
    }
    IP_one.type=basepairs_WC[ii].type;
    int index_i=Molecule::get_residue_index(IP_one.residue_i);
    int index_j=Molecule::get_residue_index(IP_one.residue_j);
    if(index_i>index_j) RNAMolecule::swap_interactionpair_type(IP_one);
    basepairs_WC[ii]=IP_one;
  }
  for(size_t ii=0;ii<basepairs_nWC.size();ii++) {
    InteractionPair IP_one;
    string res_mac_i=basepairs_nWC[ii].residue_i;
    string res_mac_j=basepairs_nWC[ii].residue_j;
    if(Mres_mca_std.count(res_mac_i))
      IP_one.residue_i=Mres_mca_std[res_mac_i];
    else {
      string str="RNA_molecule.cpp: Cannot find standard residue name for MC-Annotate resiude! --> "+res_mac_i;
      error->one(str);
    }
    if(Mres_mca_std.count(res_mac_j))
      IP_one.residue_j=Mres_mca_std[res_mac_j];
    else {
      string str="RNA_molecule.cpp: Cannot find standard residue name for MC-Annotate resiude! --> "+res_mac_j;
      error->one(str);
    }
    IP_one.type=basepairs_nWC[ii].type;
    int index_i=Molecule::get_residue_index(IP_one.residue_i);
    int index_j=Molecule::get_residue_index(IP_one.residue_j);
    if(index_i>index_j) RNAMolecule::swap_interactionpair_type(IP_one);
    basepairs_nWC[ii]=IP_one;
  }
  for(size_t ii=0;ii<basestacks.size();ii++) {
    InteractionPair IP_one;
    string res_mac_i=basestacks[ii].residue_i;
    string res_mac_j=basestacks[ii].residue_j;
    if(Mres_mca_std.count(res_mac_i))
      IP_one.residue_i=Mres_mca_std[res_mac_i];
    else {
      string str="RNA_molecule.cpp: Cannot find standard residue name for MC-Annotate resiude! --> "+res_mac_i;
      error->one(str);
    }
    if(Mres_mca_std.count(res_mac_j))
      IP_one.residue_j=Mres_mca_std[res_mac_j];
    else {
      string str="RNA_molecule.cpp: Cannot find standard residue name for MC-Annotate resiude! --> "+res_mac_j;
      error->one(str);
    }
    vector<string> stoken=split_string(IP_one.residue_i,"/");
    IP_one.type=stoken[2];
    int index_i=Molecule::get_residue_index(IP_one.residue_i);
    stoken=split_string(IP_one.residue_j,"/");
    IP_one.type+=stoken[2];
    IP_one.type+="_"+basestacks[ii].type;
    int index_j=Molecule::get_residue_index(IP_one.residue_j);
    if(index_i>index_j) RNAMolecule::swap_interactionpair_type(IP_one);
    basestacks[ii]=IP_one;
  }
}

void RNAMolecule::clear() {
  sequences.clear();
  Molecule::clear();
  basepairs_WC.clear();
  basepairs_nWC.clear();
  basestacks.clear();
}

bool RNAMolecule::is_basepaired(string nt_i, string nt_j) const {
  Residue nt_i_one=Molecule::find_residue(nt_i);
  Residue nt_j_one=Molecule::find_residue(nt_j);
  return RNAMolecule::is_basepaired(nt_i_one,nt_j_one);
}

bool RNAMolecule::is_basepaired(Residue nt_i, Residue nt_j) const {
  nt_i.recovery_RNA_backbone_atoms();
  nt_i.simplify_RNA_base();
  nt_j.recovery_RNA_backbone_atoms();
  nt_j.simplify_RNA_base();
  double dis=nt_i["C4'"].distance(nt_j["C4'"]);
  if(dis>18.0) return 0;
  vector<double> Vxyz_tar=nt_i.get_RNA_residue_xyz();
  vector<double> Vxyz=nt_j.get_RNA_residue_xyz();
  Vxyz_tar.insert(Vxyz_tar.end(),Vxyz.begin(),Vxyz.end());
  vector<Atom> bpi=init->RESIDUES_XYZ["BPi"];
  vector<Atom> bpj=init->RESIDUES_XYZ["BPj"];
  vector<double> Vxyz_ref;
  for(size_t ii=0;ii<bpi.size();ii++) {
    vector<double> xyz=bpi[ii].get_atom_xyz();
    Vxyz_ref.insert(Vxyz_ref.end(),xyz.begin(),xyz.end());
  }
  for(size_t ii=0;ii<bpj.size();ii++) {
    vector<double> xyz=bpj[ii].get_atom_xyz();
    Vxyz_ref.insert(Vxyz_ref.end(),xyz.begin(),xyz.end());
  }
  Superpose sp_one(Vxyz_tar,Vxyz_ref);
  double rmsd1=sp_one.get_RMSD();
  Vxyz_ref.clear();
  for(size_t ii=0;ii<bpj.size();ii++) {
    vector<double> xyz=bpj[ii].get_atom_xyz();
    Vxyz_ref.insert(Vxyz_ref.end(),xyz.begin(),xyz.end());
  }
  for(size_t ii=0;ii<bpi.size();ii++) {
    vector<double> xyz=bpi[ii].get_atom_xyz();
    Vxyz_ref.insert(Vxyz_ref.end(),xyz.begin(),xyz.end());
  }
  sp_one.set_x1x2(Vxyz_tar,Vxyz_ref);
  double rmsd2=sp_one.get_RMSD();

  if(rmsd1 <= init->RNA_BP_RMSD_cutoff || rmsd2 <= init->RNA_BP_RMSD_cutoff) return 1;
  else return 0;
}

void RNAMolecule::swap_interactionpair_type(InteractionPair &IP_one) {
  string residue_i=IP_one.residue_i;
  string residue_j=IP_one.residue_j;
  string type="";
  if(IP_one.type.size()==2) {
    type=to_string((char) IP_one.type[1]);
    type+=to_string((char) IP_one.type[0]);
  }
  else if(IP_one.type.size()==6) {
    type=to_string((char) IP_one.type[1]);
    type+=to_string((char) IP_one.type[0]);
    type+=to_string((char) IP_one.type[2]);
    type+=to_string((char) IP_one.type[3]);
    type+=to_string((char) IP_one.type[5]);
    type+=to_string((char) IP_one.type[4]);
  }
  else {
    if(IP_one.type.find("_")!=string::npos) {
      type=to_string((char) IP_one.type[1]);
      type+=to_string((char) IP_one.type[0]);
      type+=IP_one.type.substr(2);
    }
    else
      type=IP_one.type;
  }

  IP_one.residue_i=residue_j;
  IP_one.residue_j=residue_i;
  IP_one.type=type;
}

vector<string> RNAMolecule::read_sequence(string seq_one) {
  vector<string> Vseq;
  string seq="";

  vector<string> stoken=split_string(seq_one);
  seq_one=stoken[0];
  for(size_t ii=1;ii<stoken.size();ii++)
    seq_one+=" "+stoken[ii];
  while(seq_one.length()>0) {
    char ch=seq_one[0];
    seq_one=seq_one.substr(1);
    if(ch==' ' || ch=='-' || ch=='&') {
      if(seq.size()>0) Vseq.push_back(seq);
      seq="";
      continue;
    }
    if(ch=='A' || ch=='a') seq+="A";
    else if(ch=='G' || ch=='g') seq+="G";
    else if(ch=='C' || ch=='c') seq+="C";
    else if(ch=='U' || ch=='u' || ch=='T' || ch=='t') seq+="U";
    else if(ch=='N' || ch=='n') seq+="N";
    else if(ch=='(') { //enable (XXX) for modified nucleotide
      size_t pos=seq_one.find_first_of(")");
      if(pos==string::npos) {
        string str="RNA_molecule.cpp: Wrong RNA sequence --> "+seq_one;
        error->one(str);
      }
      seq+="("+seq_one.substr(0,pos)+")";
      seq_one=seq_one.substr(pos+1);
    }
    else {
      string str="RNA_molecule.cpp: Wrong RNA sequence (only A/a/C/c/G/g/U/u) --> "+seq_one;
      error->one(str);
    }
  }
  if(seq.size()>0) Vseq.push_back(seq);

  return Vseq;
}

vector<string> RNAMolecule::get_resnames_from_sequence(string seq) {
  vector<string> Vresname_one;
  while(seq.length()>0) {
    string resname_one;
    if(seq[0]!='(') {
      resname_one=to_string((char)seq[0]);
      seq=seq.substr(1);
      string sname="R_"+resname_one;
      if(init->SEQUENCE_MAP.count(sname)) {
        resname_one=init->SEQUENCE_MAP[sname];
      }
      Vresname_one.push_back(resname_one);
    }
    else {  //enable (XXX) for modified nucleotide
      seq=seq.substr(0);
      size_t pos=seq.find_first_of(")");
      if(pos==string::npos) {
        string str="RNA_molecule.cpp: Wrong RNA sequence --> "+seq;
        error->one(str);
      }
      resname_one=seq.substr(0,pos);
      seq=seq.substr(pos+1);
      Vresname_one.push_back(resname_one);
    }
  }

  return Vresname_one;
}

inline bool RNAMolecule::is_WC_basepair(string type) const{
  string WC[]={"AU_cWW","UA_cWW","GC_cWW","CG_cWW","GU_cWW","UG_cWW"};
  bool check;
  if(count(begin(WC),end(WC),type)==1) check=1;
  else check=0;

  return check;
}

inline string RNAMolecule::change_dssr_residue_name(string nt_one) {
  string res_one="";
  vector<string> stoken=split_string(nt_one,".");
  if(stoken.size()==3) {
    res_one=stoken[0]+"/"+stoken[2]+"/"+stoken[1];
  }
  else if(stoken.size()==2) {
    string resname_one=stoken[1].substr(0,1);
    string resid_one=stoken[1].substr(1);
    res_one=stoken[0]+"/"+resid_one+"/"+resname_one;
  }
  else {
    string str="RNA_molecule.cpp: Wrong nucleotide name  --> "+nt_one;
    error->one(str);
  }

  return res_one;
}