#include "basic.h"
#include "constant.h"
#include "initialization.h"
#include "error.h"
#include "parse_RNA_motif.h"

using namespace CBS_NS;

enum{SL,TI,MW,SLFL,TIFL,MWFL,COIL};
const int ParseRNAMotif::sl=SL;
const int ParseRNAMotif::ti=TI;
const int ParseRNAMotif::mw=MW;
const int ParseRNAMotif::slfl=SLFL;
const int ParseRNAMotif::tifl=TIFL;
const int ParseRNAMotif::mwfl=MWFL;
const int ParseRNAMotif::coil=COIL;

ParseRNAMotif::ParseRNAMotif(Cobmms *cbs) : Pointers(cbs), basepair(cbs) {}

ParseRNAMotif::ParseRNAMotif(const ParseRNAMotif &other) : Pointers(other.cbs), basepair(other.cbs) {
  motifs=other.motifs;
  basepair=other.basepair;
  helices=other.helices;
  hx_ends=other.hx_ends;
  loops=other.loops;
}

ParseRNAMotif& ParseRNAMotif::operator=(const ParseRNAMotif &other) {
  if(this != &other) {
    Pointers::operator=(other.cbs);
    motifs=other.motifs;
    basepair=other.basepair;
    helices=other.helices;
    hx_ends=other.hx_ends;
    loops=other.loops;
  }

  return *this;
}

ParseRNAMotif::~ParseRNAMotif() {
  motifs.clear();
  motifs.shrink_to_fit();
  basepair.clear();
  helices.clear();
  helices.shrink_to_fit();
  hx_ends.clear();
  loops.clear();
  loops.shrink_to_fit();
}

void ParseRNAMotif::parse_RNA_motifs(const RNABasePair &basepair_one) {
  motifs.clear();
  basepair=basepair_one;
 
  if(basepair.mapBP.size()==0) {
    string str="parse_RNA_motif.cpp: No secondary structure --> Please set the 2D structure first!";
    error->warning(str);
  }

  ParseRNAMotif::get_helices();
  if(helices.size()==0) {
    string str="parse_RNA_motif.cpp: No helix!";
    error->warning(str);
  }
  
  ParseRNAMotif::get_loops();
  if(loops.size()==0) {
    string str="parse_RNA_motif.cpp: No loop!";
    error->warning(str);
  }
  
  if(helices.size()>0) ParseRNAMotif::extract_helix();
  
  if(loops.size()==0) return;

  ParseRNAMotif::extract_tail();
  ParseRNAMotif::extract_hairpin();
  ParseRNAMotif::extract_kiss();
  ParseRNAMotif::extract_junction();
  ParseRNAMotif::extract_loop();

  for(size_t ii=0;ii<motifs.size();ii++) {
    if(motifs[ii].type.find("way_junction")==string::npos || motifs[ii].type.find("way_junction_kiss")!=string::npos)
      continue;
    for(size_t jj=ii+1;jj<motifs.size();jj++) {
      string type_j=motifs[jj].type;
      if(type_j.find("way_junction")==string::npos || type_j.find("way_junction_kiss")!=string::npos)
        continue;
      if(type_j!=motifs[ii].type)
        continue;
      RNAMotif mtf=motifs[jj];
      if(motifs[ii].is_identical_junction(mtf))
        motifs[jj]=mtf;
    }
  }
}

int ParseRNAMotif::get_RNA_topology() const {
  string topo="Stem_loop";

  string TI[]={"THelix","Tertiary_helix","Loop_kiss_loop","Hairpin_kiss","-way_junction_kiss"};
  for(size_t ii=0;ii<motifs.size();ii++) {
    string type_one=motifs[ii].type;
    for(int jj=0;jj<get_array_size(TI);jj++) {
      if(type_one.find(TI[jj])!=string::npos) {
        topo="Tertiary_interaction";
        break;
      }
    }
  }

  if(topo=="Stem_loop") {
    for(size_t ii=0;ii<motifs.size();ii++) {
      if(motifs[ii].type.find("-way_junction")!=string::npos && motifs[ii].type!="2-way_junction") {
        topo="Multi-way_junction";
        break;
      }
    }
  }

  int check_FL=0;
  for(size_t ii=0;ii<motifs.size();ii++) {
    if(motifs[ii].type=="Loop" || motifs[ii].type=="Loop_kiss_loop") {
      if(ParseRNAMotif::is_flexible_loop(motifs[ii])) {
        check_FL=1;
        break;
      }
    }
  }

  if(check_FL==1)
    topo+="_with_flexible_loop";
  
  size_t n_unpaired_chain=0;
  for(size_t ii=0;ii<motifs.size();ii++) {
    if(motifs[ii].type=="Unpaired_chain")
      n_unpaired_chain++;
  }
  if(n_unpaired_chain==motifs.size())
    topo="Coil";

  int topo_value=-1;
  if(topo=="Stem_loop") topo_value=ParseRNAMotif::sl;
  else if(topo=="Stem_loop_with_flexible_loop") topo_value=ParseRNAMotif::slfl;
  else if(topo=="Tertiary_interaction") topo_value=ParseRNAMotif::ti;
  else if(topo=="Tertiary_interaction_with_flexible_loop") topo_value=ParseRNAMotif::tifl;
  else if(topo=="Multi-way_junction") topo_value=ParseRNAMotif::mw;
  else if(topo=="Multi-way_junction_with_flexible_loop") topo_value=ParseRNAMotif::mwfl;
  else if(topo=="Coil") topo_value=ParseRNAMotif::coil;

  return topo_value;
}

string ParseRNAMotif::get_RNA_topology_string() const {
  string topo="Stem_loop";

  string TI[]={"THelix","Tertiary_helix","Loop_kiss_loop","Hairpin_kiss","-way_junction_kiss"};
  for(size_t ii=0;ii<motifs.size();ii++) {
    string type_one=motifs[ii].type;
    for(int jj=0;jj<get_array_size(TI);jj++) {
      if(type_one.find(TI[jj])!=string::npos) {
        topo="Tertiary_interaction";
        break;
      }
    }
  }

  if(topo=="Stem_loop") {
    for(size_t ii=0;ii<motifs.size();ii++) {
      if(motifs[ii].type.find("-way_junction")!=string::npos && motifs[ii].type!="2-way_junction") {
        topo="Multi-way_junction";
        break;
      }
    }
  }

  int check_FL=0;
  for(size_t ii=0;ii<motifs.size();ii++) {
    if(motifs[ii].type=="Loop" || motifs[ii].type=="Loop_kiss_loop") {
      if(ParseRNAMotif::is_flexible_loop(motifs[ii])) {
        check_FL=1;
        break;
      }
    }
  }

  if(check_FL==1)
    topo+="_with_flexible_loop";
  
  size_t n_unpaired_chain=0;
  for(size_t ii=0;ii<motifs.size();ii++) {
    if(motifs[ii].type=="Unpaired_chain")
      n_unpaired_chain++;
  }
  if(n_unpaired_chain==motifs.size())
    topo="Coil";

  return topo;
}

bool ParseRNAMotif::is_flexible_loop(const RNAMotif &mtf_one) const {
  if(mtf_one.type!="Loop" && mtf_one.type!="Loop_kiss_loop") return false;

  size_t nt_s=basepair.get_residue_index(mtf_one.ends[0]);
  size_t nt_e=basepair.get_residue_index(mtf_one.ends.back());
  string secstr_before=basepair.secstr.substr(0,nt_s+1);
  string secstr_after=basepair.secstr.substr(nt_e);

  int check_s=1;
  int nn=((int)init->RNA_BP_NOTATION.size()-1)/2;
  for(int ii=0;ii<nn;ii++) {
    char Ln=init->RNA_BP_NOTATION[ii*2+1];
    char Rn=init->RNA_BP_NOTATION[ii*2+2];
    int nL=0,nR=0;
    for(size_t jj=0;jj<secstr_before.size();jj++) {
      if(secstr_before[jj]==Ln) nL++;
      else if(secstr_before[jj]==Rn) nR++;
    }
    if(nL!=nR) {
      check_s=0;
      break;
    }
  }
  int check_e=1;
  for(int ii=0;ii<nn;ii++) {
    char Ln=init->RNA_BP_NOTATION[ii*2+1];
    char Rn=init->RNA_BP_NOTATION[ii*2+2];
    int nL=0,nR=0;
    for(size_t jj=0;jj<secstr_after.size();jj++) {
      if(secstr_after[jj]==Ln) nL++;
      else if(secstr_after[jj]==Rn) nR++;
    }
    if(nL!=nR) {
      check_e=0;
      break;
    }
  }

  if(check_s==1 || check_e==1) return true;
  else return false;
}

size_t ParseRNAMotif::count_flexible_loop() const {
  size_t num=0;
  for(size_t ii=0;ii<motifs.size();ii++) {
    if(motifs[ii].type=="Loop" || motifs[ii].type=="Loop_kiss_loop")
      if(ParseRNAMotif::is_flexible_loop(motifs[ii])) num++;
  }

  return num;
}

size_t ParseRNAMotif::count_unpaired_chain() const {
  size_t num=0;
  for(size_t ii=0;ii<motifs.size();ii++)
    if(motifs[ii].type=="Unpaired_chain") num++;

  return num;
}

size_t ParseRNAMotif::size() const {
  return motifs.size();
}

void ParseRNAMotif::clear() {
  motifs.clear();
  motifs.shrink_to_fit();
  basepair.clear();
  helices.clear();
  helices.shrink_to_fit();
  hx_ends.clear();
  loops.clear();
  loops.shrink_to_fit();
}

void ParseRNAMotif::get_helices() {
  helices.clear();
  hx_ends.clear();
  map<string,string> BPs=basepair.mapBP;
  while(BPs.size()>0) {
    Chain hx_c1(cbs),hx_c2(cbs);
    string bpi=BPs.begin()->first;
    string bpj=BPs.begin()->second;
    BPs.erase(bpi);
    hx_c1.insert_residue(bpi,Chain::begin);
    hx_c2.insert_residue(bpj,Chain::end);
    int bp_level=basepair.get_basepair_level(bpi,bpj);
    string bpi_previous=bpi;
    string bpj_next=bpj;
    while (1) {
      bpi_previous=basepair.previous_residue(bpi_previous);
      bpj_next=basepair.next_residue(bpj_next);
      //cout<<"i/j/i-1/j+1: "<<bpi<<" "<<bpj<<" "<<bpi_previous<<" "<<bpj_next<<endl;
      if(bpi_previous==ResidueConst::RESIDUE_NULL || bpj_next==ResidueConst::RESIDUE_NULL) break;
      if(BPs.count(bpi_previous)==0 && BPs.count(bpj_next)==0) break;
      if(basepair.is_basepaired(bpi_previous,bpj_next)) {
        hx_c1.insert_residue(bpi_previous,Chain::begin);
        hx_c2.insert_residue(bpj_next,Chain::end);
        if(BPs.count(bpi_previous)) BPs.erase(bpi_previous);
        if(BPs.count(bpj_next)) BPs.erase(bpj_next);
      }
      else if(bp_level>1) { //enable 1 unpaired nucleotide in each side of tertiary helix
        string bpi_previous_1=basepair.previous_residue(bpi_previous);
        string bpj_next_1=basepair.next_residue(bpj_next);
        if(bpj_next_1!=ResidueConst::RESIDUE_NULL && BPs.count(bpi_previous) && basepair.is_basepaired(bpi_previous,bpj_next_1)) {
          hx_c1.insert_residue(bpi_previous,Chain::begin);
          hx_c2.insert_residue(bpj_next,Chain::end);
          hx_c2.insert_residue(bpj_next_1,Chain::end);
          bpj_next=bpj_next_1;
          if(BPs.count(bpi_previous)) BPs.erase(bpi_previous);
          if(BPs.count(bpj_next)) BPs.erase(bpj_next);
        }
        else if(bpi_previous_1!=ResidueConst::RESIDUE_NULL && BPs.count(bpi_previous_1) && basepair.is_basepaired(bpi_previous_1,bpj_next)) {
          hx_c1.insert_residue(bpi_previous,Chain::begin);
          hx_c1.insert_residue(bpi_previous_1,Chain::begin);
          hx_c2.insert_residue(bpj_next,Chain::end);
          bpi_previous=bpi_previous_1;
          if(BPs.count(bpi_previous)) BPs.erase(bpi_previous);
          if(BPs.count(bpj_next)) BPs.erase(bpj_next);
        }
        else if(bpi_previous_1!=ResidueConst::RESIDUE_NULL && BPs.count(bpi_previous_1) && bpj_next_1!=ResidueConst::RESIDUE_NULL && basepair.is_basepaired(bpi_previous_1,bpj_next_1)) {
          hx_c1.insert_residue(bpi_previous,Chain::begin);
          hx_c1.insert_residue(bpi_previous_1,Chain::begin);
          hx_c2.insert_residue(bpj_next,Chain::end);
          hx_c2.insert_residue(bpj_next_1,Chain::end);
          bpi_previous=bpi_previous_1;
          bpj_next=bpj_next_1;
          if(BPs.count(bpi_previous)) BPs.erase(bpi_previous);
          if(BPs.count(bpj_next)) BPs.erase(bpj_next);
        }
        else break;
      }
      else break;
    }

    string bpi_next=bpi;
    string bpj_previous=bpj;
    while (1) {
      bpi_next=basepair.next_residue(bpi_next);
      bpj_previous=basepair.previous_residue(bpj_previous);
      //cout<<"i/j/i+1/j-1: "<<bpi<<" "<<bpj<<" "<<bpi_next<<" "<<bpj_previous<<endl;
      if(bpi_next==ResidueConst::RESIDUE_NULL || bpj_previous==ResidueConst::RESIDUE_NULL) break;
      if(BPs.count(bpi_next)==0 && BPs.count(bpj_previous)==0) break;
      if(basepair.is_basepaired(bpi_next,bpj_previous)) {
        hx_c1.insert_residue(bpi_next,Chain::end);
        hx_c2.insert_residue(bpj_previous,Chain::begin);
        if(BPs.count(bpi_next)) BPs.erase(bpi_next);
        if(BPs.count(bpj_previous)) BPs.erase(bpj_previous);
      }
      else if(bp_level>1) { //enable 1 unpaired nucleotide in each side of tertiary helix
        string bpi_next_1=basepair.next_residue(bpi_next);
        string bpj_previous_1=basepair.previous_residue(bpj_previous);
        if(bpj_previous_1!=ResidueConst::RESIDUE_NULL && BPs.count(bpi_next) && basepair.is_basepaired(bpi_next,bpj_previous_1)) {
          hx_c1.insert_residue(bpi_next,Chain::end);
          hx_c2.insert_residue(bpj_previous,Chain::begin);
          hx_c2.insert_residue(bpj_previous_1,Chain::begin);
          bpj_previous=bpj_previous_1;
          if(BPs.count(bpi_next)) BPs.erase(bpi_next);
          if(BPs.count(bpj_previous)) BPs.erase(bpj_previous);
        }
        else if(bpi_next_1!=ResidueConst::RESIDUE_NULL && BPs.count(bpi_next_1) && basepair.is_basepaired(bpi_next_1,bpj_previous)) {
          hx_c1.insert_residue(bpi_next,Chain::end);
          hx_c1.insert_residue(bpi_next_1,Chain::end);
          hx_c2.insert_residue(bpj_previous,Chain::begin);
          bpi_next=bpi_next_1;
          if(BPs.count(bpi_next)) BPs.erase(bpi_next);
          if(BPs.count(bpj_previous)) BPs.erase(bpj_previous);
        }
        else if(bpi_next_1!=ResidueConst::RESIDUE_NULL && BPs.count(bpi_next_1) && bpj_previous_1!=ResidueConst::RESIDUE_NULL && basepair.is_basepaired(bpi_next_1,bpj_previous_1)) {
          hx_c1.insert_residue(bpi_next,Chain::end);
          hx_c1.insert_residue(bpi_next_1,Chain::end);
          hx_c2.insert_residue(bpj_previous,Chain::begin);
          hx_c2.insert_residue(bpj_previous_1,Chain::begin);
          bpi_next=bpi_next_1;
          bpj_previous=bpj_previous_1;
          if(BPs.count(bpi_next)) BPs.erase(bpi_next);
          if(BPs.count(bpj_previous)) BPs.erase(bpj_previous);
        }
        else break;
      }
      else break;
    }

    if(hx_c1.size()>1 && hx_c2.size()>1) {
      Molecule hx(cbs);
      hx_c1.chain_label=Chain::rna;
      hx_c2.chain_label=Chain::rna;
      hx.insert_chain(hx_c1);
      hx.insert_chain(hx_c2);
      helices.push_back(hx);
      string bp_s=hx_c1.first_of_residue();
      string bp_e=hx_c2.last_of_residue();
      hx_ends.insert(pair<string,string>(bp_s,bp_e));
      bp_s=hx_c2.first_of_residue();
      bp_e=hx_c1.last_of_residue();
      hx_ends.insert(pair<string,string>(bp_s,bp_e));
    }
  }
}

void ParseRNAMotif::get_loops() {
  loops.clear();
  for(auto& it : hx_ends) {
    string res_i=it.second;
    Chain loop_c(cbs);
    loop_c.insert_residue(res_i,Chain::begin);
    while(1) {
      res_i=basepair.next_residue(res_i);
      if(res_i==ResidueConst::RESIDUE_NULL) break;
      loop_c.insert_residue(res_i,Chain::end);
      if(hx_ends.count(res_i)) break;
      if(res_i==it.second) break;
    }
    if(loop_c.size()>1) {
      Molecule loop(cbs);
      loop_c.chain_label=Chain::rna;
      loop.insert_chain(loop_c);
      loops.push_back(loop);
    }
  }

  vector<string> unsigned_res;
  for(size_t ii=0;ii<basepair.Vresnames.size();ii++) {
    string res_i=basepair.Vresnames[ii];
    int check=0;
    for(size_t jj=0;jj<helices.size();jj++) {
      if(count(helices[jj].Vresnames.begin(),helices[jj].Vresnames.end(),res_i)) {
        check=1;
        break;
      }
    }
    if(!check) {
      for(size_t kk=0;kk<loops.size();kk++) {
        if(count(loops[kk].Vresnames.begin(),loops[kk].Vresnames.end(),res_i)) {
          check=1;
          break;
        }
      }
    }
    if(!check) unsigned_res.push_back(res_i);
  }

  while(unsigned_res.size()>0) {
    Chain loop_c(cbs);
    string res_i=unsigned_res[0];
    loop_c.insert_residue(res_i,Chain::begin);
    unsigned_res.erase(unsigned_res.begin());
    while(1) {
      res_i=basepair.next_residue(res_i);
      if(res_i==ResidueConst::RESIDUE_NULL) break;
      if(res_i==loop_c.first_of_residue()) break;
      loop_c.insert_residue(res_i,Chain::end);
      if(hx_ends.count(res_i)) break;
      for(size_t jj=0;jj<unsigned_res.size();jj++) {
        if(res_i==unsigned_res[jj]) {
          unsigned_res.erase(unsigned_res.begin()+jj);
          break;
        }
      }
    }
    Molecule loop(cbs);
    loop_c.chain_label=Chain::rna;
    loop.insert_chain(loop_c);
    loops.push_back(loop);
  }

}

void ParseRNAMotif::extract_helix() {
  vector<Molecule> Thelix=helices;
  helices.clear();
  for(size_t ii=0;ii<Thelix.size();ii++) {
    string bp_s=Thelix[ii].chains[0].first_of_residue();
    string bp_e=Thelix[ii].chains[1].last_of_residue();
    int bp_level=basepair.get_basepair_level(bp_s,bp_e);
    if(bp_level==1) {
      RNAMotif hx(cbs);
      hx.type="Helix";
      hx.insert_RNA_segment(Thelix[ii].chains[0],0);
      string secstr_one=basepair.get_segment_secondary_structure(bp_s,Thelix[ii].chains[0].last_of_residue());
      hx.secstrs.push_back(secstr_one);
      hx.insert_RNA_segment(Thelix[ii].chains[1],0);
      secstr_one=basepair.get_segment_secondary_structure(Thelix[ii].chains[1].first_of_residue(),bp_e);
      hx.secstrs.push_back(secstr_one);
      motifs.push_back(hx);
    }
    else {
      helices.push_back(Thelix[ii]);
    }
  }
  Thelix.clear();

  for(size_t ii=0;ii<helices.size();ii++) {
    RNAMotif thx(cbs);
    if(helices[ii].chains[0].size()==helices[ii].chains[1].size()) {
      thx.type="THelix";
      thx.insert_RNA_segment(helices[ii].chains[0],0);
      string secstr_one=basepair.get_segment_secondary_structure(helices[ii].chains[0].first_of_residue(),helices[ii].chains[0].last_of_residue());
      thx.secstrs.push_back(secstr_one);
      thx.insert_RNA_segment(helices[ii].chains[1],0);
      secstr_one=basepair.get_segment_secondary_structure(helices[ii].chains[1].first_of_residue(),helices[ii].chains[1].last_of_residue());
      thx.secstrs.push_back(secstr_one);
      motifs.push_back(thx);
    }
    else {
      Chain hx_c1(cbs), hx_c2(cbs);
      vector<Chain> Vhx_c1,Vhx_c2;
      int index_c1=0;
      int index_c2=static_cast<int>(helices[ii].chains[1].size()-1);
      while(index_c1<(int)helices[ii].chains[0].size() && index_c2>=0) {
        string bp_i=helices[ii].chains[0].residues[index_c1].get_residue_name();
        string bp_j=helices[ii].chains[1].residues[index_c2].get_residue_name();
        if(basepair.is_basepaired(bp_i,bp_j)) {
          hx_c1.insert_residue(bp_i,Chain::end);
          hx_c2.insert_residue(bp_j,Chain::begin);
          index_c1++;
          index_c2--;
        }
        else if(basepair.is_basepaired(bp_i) && !basepair.is_basepaired(bp_j)) {
          index_c2--;
          if(hx_c1.size()>1) {
            hx_c1.chain_label=Chain::rna;
            hx_c2.chain_label=Chain::rna;
            Vhx_c1.push_back(hx_c1);
            Vhx_c2.push_back(hx_c2);
          }
          hx_c1.clear();
          hx_c2.clear();
        }
        else if(!basepair.is_basepaired(bp_i) && basepair.is_basepaired(bp_j)) {
          index_c1++;
          if(hx_c1.size()>1) {
            hx_c1.chain_label=Chain::rna;
            hx_c2.chain_label=Chain::rna;
            Vhx_c1.push_back(hx_c1);
            Vhx_c2.push_back(hx_c2);
          }
          hx_c1.clear();
          hx_c2.clear();
        }
        else {
          index_c1++;
          index_c2--;
          if(hx_c1.size()>1) {
            hx_c1.chain_label=Chain::rna;
            hx_c2.chain_label=Chain::rna;
            Vhx_c1.push_back(hx_c1);
            Vhx_c2.push_back(hx_c2);
          }
          hx_c1.clear();
          hx_c2.clear();
        }
      }

      for(size_t jj=0;jj<Vhx_c1.size();jj++) {
        thx.type="THelix";
        thx.insert_RNA_segment(Vhx_c1[jj],0);
        string secstr_one=basepair.get_segment_secondary_structure(Vhx_c1[jj].first_of_residue(),Vhx_c1[jj].last_of_residue());
        thx.secstrs.push_back(secstr_one);
        thx.insert_RNA_segment(Vhx_c2[jj],0);
        secstr_one=basepair.get_segment_secondary_structure(Vhx_c2[jj].first_of_residue(),Vhx_c2[jj].last_of_residue());
        thx.secstrs.push_back(secstr_one);
        motifs.push_back(thx);
        thx.clear();
      }
    }
  }
}

void ParseRNAMotif::extract_tail() {
  vector<Molecule> Tloop=loops;
  loops.clear();
  for(size_t ii=0;ii<Tloop.size();ii++) {
    string res_s=Tloop[ii].chains[0].first_of_residue();
    string res_e=Tloop[ii].chains[0].last_of_residue();
    string res_s_p=basepair.previous_residue(res_s);
    string res_e_n=basepair.next_residue(res_e);
    string secstr_one=basepair.get_segment_secondary_structure(res_s,res_e);
    RNAMotif mtf(cbs);
    if((res_s_p==ResidueConst::RESIDUE_NULL && res_e_n==ResidueConst::RESIDUE_NULL) || (res_s_p==res_e && res_e_n==res_s)) {
      mtf.type="Unpaired_chain";
      mtf.insert_RNA_segment(Tloop[ii].chains[0],0);
      mtf.secstrs.push_back(secstr_one);
      motifs.push_back(mtf);
    }
    else if(res_s_p==ResidueConst::RESIDUE_NULL && res_e_n!=ResidueConst::RESIDUE_NULL) {
      mtf.type="Tail5";
      mtf.insert_RNA_segment(Tloop[ii].chains[0],-1);
      mtf.secstrs.push_back(secstr_one);
      motifs.push_back(mtf);
    }
    else if(res_s_p!=ResidueConst::RESIDUE_NULL && res_e_n==ResidueConst::RESIDUE_NULL) {
      mtf.type="Tail3";
      mtf.insert_RNA_segment(Tloop[ii].chains[0],-1);
      mtf.secstrs.push_back(secstr_one);
      motifs.push_back(mtf);
    }
    else
      loops.push_back(Tloop[ii]);
  }
}

void ParseRNAMotif::extract_hairpin() {
  vector<Molecule> Tloop=loops;
  loops.clear();
  for(size_t ii=0;ii<Tloop.size();ii++) {
    string res_s=Tloop[ii].chains[0].first_of_residue();
    string res_e=Tloop[ii].chains[0].last_of_residue();
    string secstr_one=basepair.get_segment_secondary_structure(res_s,res_e);
    if(basepair.is_basepaired(res_s,res_e)) {
      RNAMotif hp(cbs);
      hp.type="Hairpin";
      hp.insert_RNA_segment(Tloop[ii].chains[0],-2);
      hp.secstrs.push_back(secstr_one);
      motifs.push_back(hp);
    }
    else
      loops.push_back(Tloop[ii]);
  }
}

void ParseRNAMotif::extract_kiss() {
  for(size_t ii=0;ii<helices.size();ii++) {
    for(int jj=0;jj<2;jj++) {
      string res_s=helices[ii].chains[jj].first_of_residue();
      string res_e=helices[ii].chains[jj].last_of_residue();
      int bp_level=basepair.get_basepair_level(res_s,hx_ends[res_s]);
      if(bp_level==1) continue;
      hx_ends.erase(res_s);

      Chain loop_s(cbs),loop_e(cbs);
      int check_s=0;
      int check_e=0;
      for(size_t kk=0;kk<loops.size();kk++) {
        string res_e_l=loops[kk].chains[0].last_of_residue();
        if(res_e_l==res_s) {
          check_s=1;
          loop_s.set_chain(loops[kk].chains[0]);
          loops.erase(loops.begin()+kk);
        }
      }
      for(size_t kk=0;kk<loops.size();kk++) {
        string res_s_l=loops[kk].chains[0].first_of_residue();
        if(res_s_l==res_e) {
          check_e=1;
          loop_e.set_chain(loops[kk].chains[0]);
          loops.erase(loops.begin()+kk);
        }
      }

      if(check_s && check_e) {
        RNAMotif mtf(cbs);
        string res_s_m=loop_s.first_of_residue();
        string res_e_m=loop_e.last_of_residue();
        if(basepair.is_basepaired(res_s_m,res_e_m)) mtf.type="Hairpin_kiss";
        else mtf.type="Loop_kiss_loop";
        Chain loop_f(cbs);
        loop_f.set_chain(loop_s);
        for(size_t ll=1;ll<helices[ii].chains[jj].size()-1;ll++)
          loop_f.insert_residue(helices[ii].chains[jj].residues[ll],Chain::end);
        for(size_t ll=0;ll<loop_e.residues.size();ll++)
          loop_f.insert_residue(loop_e.residues[ll],Chain::end);
        mtf.insert_RNA_segment(loop_f,-2);
        mtf.lengths.pop_back();
        int len=(int)loop_s.size()-2;
        mtf.lengths.push_back(len);
        len=(int)helices[ii].chains[jj].size();
        mtf.lengths.push_back(len);
        len=(int)loop_e.size()-2;
        mtf.lengths.push_back(len);
        string res_s_f=loop_f.first_of_residue();
        string res_e_f=loop_f.last_of_residue();
        string secstr_one=basepair.get_segment_secondary_structure(res_s_f,res_e_f);
        mtf.secstrs.push_back(secstr_one);
        motifs.push_back(mtf);
      }
      else {
        RNAMotif mtf(cbs);
        mtf.type="Tertiary_helix";
        mtf.insert_RNA_segment(helices[ii].chains[jj],0);
        string res_s_h=helices[ii].chains[jj].first_of_residue();
        string res_e_h=helices[ii].chains[jj].last_of_residue();
        string secstr_one=basepair.get_segment_secondary_structure(res_s_h,res_e_h);
        mtf.secstrs.push_back(secstr_one);
        motifs.push_back(mtf);
        if(check_s) {
          RNAMotif lp(cbs);
          lp.type="Loop";
          lp.insert_RNA_segment(loop_s,-2);
          string res_s_l=loop_s.first_of_residue();
          string res_e_l=loop_s.last_of_residue();
          string secstr_one=basepair.get_segment_secondary_structure(res_s_l,res_e_l);
          lp.secstrs.push_back(secstr_one);
          motifs.push_back(lp);
        }
        else if(check_e) {
          RNAMotif lp(cbs);
          lp.type="Loop";
          lp.insert_RNA_segment(loop_e,-2);
          string res_s_l=loop_e.first_of_residue();
          string res_e_l=loop_e.last_of_residue();
          string secstr_one=basepair.get_segment_secondary_structure(res_s_l,res_e_l);
          lp.secstrs.push_back(secstr_one);
          motifs.push_back(lp);
        }
      }
    }
  }
}

void ParseRNAMotif::extract_junction() {
  for(auto& it : hx_ends) {
    map<string,int> hx_end_loop;
    for(auto& it1 : hx_ends) {
      int nloop=-1;
      for(size_t ii=0;ii<loops.size();ii++) {
        string res_s=loops[ii].chains[0].first_of_residue();
        if(res_s==it1.second) {
          nloop=ii;
          break;
        }
      }
      if(nloop>=0) hx_end_loop.insert(pair<string,int>(it1.second,nloop));

      nloop=-1;
      for(size_t ii=0;ii<loops.size();ii++) {
        string res_e=loops[ii].chains[0].last_of_residue();
        if(res_e==it1.first) {
          nloop=ii;
          break;
        }
      }
      if(nloop>=0) hx_end_loop.insert(pair<string,int>(it1.first,nloop));
    }

    vector<int> Vloop;
    if(hx_end_loop.count(it.second)==0) continue;
    Vloop.push_back(hx_end_loop[it.second]);
    hx_end_loop.erase(it.second);
    while(1) {
      int check1=0;
      for(auto& it1 : hx_ends) {
        if(it1.first==it.first) continue;
        if(hx_end_loop.count(it1.first)==0 || hx_end_loop.count(it1.second)==0) continue;
        if(hx_end_loop[it1.first]==Vloop.back()) {
          Vloop.push_back(hx_end_loop[it1.second]);
          check1=1;
          hx_end_loop.erase(it1.first);
          hx_end_loop.erase(it1.second);
          break;
        }
      }
      if(check1==0) break;
    }

    if(hx_end_loop.count(it.first) && Vloop.back()==hx_end_loop[it.first]) {
      size_t ns=0;
      int nmin=loops[Vloop[0]].chains[0].size();
      for(size_t jj=1;jj<Vloop.size();jj++) {
        int len=loops[Vloop[jj]].chains[0].size();
        if(len<nmin) {
          nmin=len;
          ns=jj;
        }
      }
      RNAMotif mtf(cbs);
      mtf.type=to_string((int)Vloop.size())+"-way_junction";
      for(size_t jj=ns;jj<Vloop.size();jj++) {
        mtf.insert_RNA_segment(loops[Vloop[jj]].chains[0],-2);
        string res_s_l=loops[Vloop[jj]].chains[0].first_of_residue();
        string res_e_l=loops[Vloop[jj]].chains[0].last_of_residue();
        string secstr_one=basepair.get_segment_secondary_structure(res_s_l,res_e_l);
        mtf.secstrs.push_back(secstr_one);
      }
      for(size_t jj=0;jj<ns;jj++) {
        mtf.insert_RNA_segment(loops[Vloop[jj]].chains[0],-2);
        string res_s_l=loops[Vloop[jj]].chains[0].first_of_residue();
        string res_e_l=loops[Vloop[jj]].chains[0].last_of_residue();
        string secstr_one=basepair.get_segment_secondary_structure(res_s_l,res_e_l);
        mtf.secstrs.push_back(secstr_one);
      }
      motifs.push_back(mtf);

      vector<Molecule> Tloop=loops;
      loops.clear();
      for(size_t jj=0;jj<Tloop.size();jj++) {
        if(count(Vloop.begin(),Vloop.end(),jj)) continue;
        else loops.push_back(Tloop[jj]);
      }
      hx_end_loop.erase(it.first);
    }
    else {
      string res_s_m=loops[Vloop[0]].chains[0].first_of_residue();
      string res_e_m=loops[Vloop.back()].chains[0].last_of_residue();
      size_t index_jj;
      int flag_check=0;
      for(size_t jj=0;jj<motifs.size();jj++) {
        if(motifs[jj].type!="Loop_kiss_loop") continue;
        string res_s_k=motifs[jj].chains[0].first_of_residue();
        string res_e_k=motifs[jj].chains[0].last_of_residue();
        if(basepair.is_basepaired(res_s_k,res_e_m) && basepair.is_basepaired(res_e_k,res_s_m)) {
          RNAMotif mtf(cbs);
          mtf.type=to_string((int)Vloop.size()+1)+"-way_junction_kiss";
          for(size_t kk=0;kk<Vloop.size();kk++) {
            mtf.insert_RNA_segment(loops[Vloop[kk]].chains[0],-2);
            string res_s_l=loops[Vloop[kk]].chains[0].first_of_residue();
            string res_e_l=loops[Vloop[kk]].chains[0].last_of_residue();
            string secstr_one=basepair.get_segment_secondary_structure(res_s_l,res_e_l);
            mtf.secstrs.push_back(secstr_one);
          }
          //int len_LKL=motifs[jj].lengths[0]+motifs[jj].lengths[1]+motifs[jj].lengths[2];
          //mtf.lengths.push_back(len_LKL);
          //mtf.ends.push_back(motifs[jj].ends[0]);
          //mtf.ends.push_back(motifs[jj].ends.back());
          //mtf.sequences.push_back(motifs[jj].sequences[0]);
          mtf.insert_RNA_segment(motifs[jj].chains[0],-2);
          mtf.secstrs.push_back(motifs[jj].secstrs[0]);
          motifs.push_back(mtf);
          index_jj=jj;
          flag_check=1;

          vector<Molecule> Tloop=loops;
          loops.clear();
          for(size_t jj=0;jj<Tloop.size();jj++) {
            if(count(Vloop.begin(),Vloop.end(),jj)) continue;
            else loops.push_back(Tloop[jj]);
          }
        }
      }
      if(flag_check)
        motifs.erase(motifs.begin()+index_jj);
    }
  }
}

void ParseRNAMotif::extract_loop() {
  for(size_t ii=0;ii<loops.size();ii++) {
    RNAMotif lp(cbs);
    lp.type="Loop";
    lp.insert_RNA_segment(loops[ii].chains[0],-2);
    string res_s_l=loops[ii].chains[0].first_of_residue();
    string res_e_l=loops[ii].chains[0].last_of_residue();
    string secstr_one=basepair.get_segment_secondary_structure(res_s_l,res_e_l);
    lp.secstrs.push_back(secstr_one);
    motifs.push_back(lp);
  }
}

ostream& operator <<(ostream &output, ParseRNAMotif &PRmotif_one) {
  output<<"List of "<<PRmotif_one.size()<<" motifs:"<<endl;
  for(size_t ii=0;ii<PRmotif_one.size();ii++) {
    string mtf_name=PRmotif_one.motifs[ii].get_RNA_motif_info();
    output<<mtf_name<<endl;
  }

  return output;
}