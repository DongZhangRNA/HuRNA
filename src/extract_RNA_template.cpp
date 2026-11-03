#include "basic.h"
#include "constant.h"
#include "readfile.h"
#include "initialization.h"
#include "error.h"
#include "RNA_molecule.h"
#include "superpose.h"
#include "extract_RNA_template.h"

using namespace CBS_NS;

const int Nw_extend_max=8;

ExtractRNATemplate::ExtractRNATemplate(Cobmms *cbs) : Pointers(cbs) {
  basepair = nullptr;
}

ExtractRNATemplate::ExtractRNATemplate(const ExtractRNATemplate &other) : Pointers(other.cbs) {
  ExtractRNATemplate::clear();
  templates=other.templates;
  extended_templates=other.extended_templates;
  Vinfo=other.Vinfo;
  helices=other.helices;
  hx_ends=other.hx_ends;
  pseudo_basepairs=other.pseudo_basepairs;
  loops=other.loops;

  basepair = nullptr;
}

ExtractRNATemplate& ExtractRNATemplate::operator=(const ExtractRNATemplate &other) {
  if(this != &other) {
    Pointers::operator=(other.cbs);
    templates=other.templates;
    extended_templates=other.extended_templates;
    Vinfo=other.Vinfo;
    helices=other.helices;
    hx_ends=other.hx_ends;
    pseudo_basepairs=other.pseudo_basepairs;
    loops=other.loops;

    basepair = nullptr;
  }

  return *this;
}

ExtractRNATemplate::~ExtractRNATemplate() {
  templates.clear();
  templates.shrink_to_fit();
  extended_templates.clear();
  extended_templates.shrink_to_fit();
  Vinfo.clear();
  Vinfo.shrink_to_fit();
  helices.clear();
  helices.shrink_to_fit();
  hx_ends.clear();
  pseudo_basepairs.clear();
  loops.clear();
  loops.shrink_to_fit();

  basepair = nullptr;
}

//flag_extend=0, no extened templates; flag_extend=1, extend templates by +/-1 closing basepair;
void ExtractRNATemplate::extract_RNA_templates(const RNAMolecule &mol_one, const RNABasePair &basepair_one, int flag_extend) {
  templates.clear();
  extended_templates.clear();
  Vinfo.clear();

  if(mol_one.Vresnames.size()==0) {
    string str="extract_RNA_template.cpp: Please set RNA molecule first!";
    error->one(str);
  }

  basepair = new RNABasePair(cbs);
  basepair->set_RNA_basepairs(basepair_one);

  if(basepair->size()==0) {
    string str="extract_RNA_template.cpp: No secondary structure!";
    Vinfo.push_back(str);
  }

  ExtractRNATemplate::get_helices();
  if(helices.size()==0) {
    string str="extract_RNA_template.cpp: No helix!";
    Vinfo.push_back(str);
    return;
  }
  ExtractRNATemplate::get_loops();
  if(loops.size()==0) {
    string str="extract_RNA_template.cpp: No loop!";
    Vinfo.push_back(str);
    return;
  }

  if(helices.size()>0) ExtractRNATemplate::extract_helix();

  pseudo_basepairs.clear();
  if(flag_extend) {
    for(auto& it : hx_ends) {
      string nt_previous=mol_one.previous_residue(it.first);
      string nt_next=mol_one.next_residue(it.second);
      if(nt_previous==ResidueConst::RESIDUE_NULL || nt_next==ResidueConst::RESIDUE_NULL) continue;
      if(mol_one.is_basepaired(nt_previous,nt_next)) {
        if(pseudo_basepairs.count(nt_previous)==0)
          pseudo_basepairs.emplace(nt_previous,nt_next);
        if(pseudo_basepairs.count(nt_next)==0)
          pseudo_basepairs.emplace(nt_next,nt_previous);
      }
    }
  }

  ExtractRNATemplate::extract_tail();
  ExtractRNATemplate::extract_hairpin();
  ExtractRNATemplate::extract_loop_template();
  ExtractRNATemplate::extract_pseudoknot_template();
  ExtractRNATemplate::extract_kiss_template();
  ExtractRNATemplate::extract_junction_template();

  extended_templates.clear();
  if(flag_extend) {
    for(size_t ii=0;ii<templates.size();ii++) {
      string type=templates[ii].type;
      if(type=="Hairpin_kiss") ExtractRNATemplate::extend_hairpin_kiss(templates[ii]);
      else if(type=="Loop_kiss_loop") ExtractRNATemplate::extend_loop_kiss_loop(templates[ii]);
      else if(type.find("-way_junction")!=string::npos) ExtractRNATemplate::extend_junction(templates[ii]);
    }
  }

  delete basepair;
  basepair=nullptr;
}

size_t ExtractRNATemplate::size() const {
  return templates.size();
}

void ExtractRNATemplate::clear() {
  templates.clear();
  templates.shrink_to_fit();
  extended_templates.clear();
  extended_templates.shrink_to_fit();
  Vinfo.clear();
  Vinfo.shrink_to_fit();
  helices.clear();
  helices.shrink_to_fit();
  hx_ends.clear();
  pseudo_basepairs.clear();
  loops.clear();
  loops.shrink_to_fit();

  basepair = nullptr;
}

void ExtractRNATemplate::get_helices() {
  helices.clear();
  hx_ends.clear();
  map<string,string> BPs=basepair->mapBP;
  while(BPs.size()>0) {
    Chain hx_c1(cbs),hx_c2(cbs);
    string bpi=BPs.begin()->first;
    string bpj=BPs.begin()->second;
    BPs.erase(bpi);
    hx_c1.insert_residue(bpi,Chain::begin);
    hx_c2.insert_residue(bpj,Chain::end);
    //cout<<"get_helices: "<<bpi<<" "<<bpj<<endl;
    int bp_level=basepair->get_basepair_level(bpi,bpj);
    string bpi_previous=bpi;
    string bpj_next=bpj;
    while (1) {
      bpi_previous=basepair->previous_residue(bpi_previous);
      bpj_next=basepair->next_residue(bpj_next);
      if(bpi_previous==ResidueConst::RESIDUE_NULL || bpj_next==ResidueConst::RESIDUE_NULL) break;
      if(BPs.count(bpi_previous)==0) break;
      if(basepair->is_basepaired(bpi_previous,bpj_next)) {
        hx_c1.insert_residue(bpi_previous,Chain::begin);
        hx_c2.insert_residue(bpj_next,Chain::end);
        BPs.erase(bpi_previous);
      }
      else if(bp_level>1) { //enable 1 unpaired nucleotide in each side of tertiary helix
        string bpi_previous_1=basepair->previous_residue(bpi_previous);
        string bpj_next_1=basepair->next_residue(bpj_next);
        if(bpj_next_1!=ResidueConst::RESIDUE_NULL && BPs.count(bpi_previous) && basepair->is_basepaired(bpi_previous,bpj_next_1)) {
          hx_c1.insert_residue(bpi_previous,Chain::begin);
          hx_c2.insert_residue(bpj_next,Chain::end);
          hx_c2.insert_residue(bpj_next_1,Chain::end);
          bpj_next=bpj_next_1;
          BPs.erase(bpi_previous);
        }
        else if(bpi_previous_1!=ResidueConst::RESIDUE_NULL && BPs.count(bpi_previous_1) && basepair->is_basepaired(bpi_previous_1,bpj_next)) {
          hx_c1.insert_residue(bpi_previous,Chain::begin);
          hx_c1.insert_residue(bpi_previous_1,Chain::begin);
          hx_c2.insert_residue(bpj_next,Chain::end);
          bpi_previous=bpi_previous_1;
          BPs.erase(bpi_previous);
        }
        else if(bpi_previous_1!=ResidueConst::RESIDUE_NULL && bpj_next_1!=ResidueConst::RESIDUE_NULL && BPs.count(bpi_previous_1) && basepair->is_basepaired(bpi_previous_1,bpj_next_1)) {
          hx_c1.insert_residue(bpi_previous,Chain::begin);
          hx_c1.insert_residue(bpi_previous_1,Chain::begin);
          hx_c2.insert_residue(bpj_next,Chain::end);
          hx_c2.insert_residue(bpj_next_1,Chain::end);
          bpi_previous=bpi_previous_1;
          bpj_next=bpj_next_1;
          BPs.erase(bpi_previous);
        }
        else break;
      }
      else break;
    }

    string bpi_next=bpi;
    string bpj_previous=bpj;
    while (1) {
      bpi_next=basepair->next_residue(bpi_next);
      bpj_previous=basepair->previous_residue(bpj_previous);
      if(bpi_next==ResidueConst::RESIDUE_NULL || bpj_previous==ResidueConst::RESIDUE_NULL) break;
      if(BPs.count(bpi_next)==0) break;
      if(basepair->is_basepaired(bpi_next,bpj_previous)) {
        hx_c1.insert_residue(bpi_next,Chain::end);
        hx_c2.insert_residue(bpj_previous,Chain::begin);
        BPs.erase(bpi_next);
      }
      else if(bp_level>1) { //enable 1 unpaired nucleotide in each side of tertiary helix
        string bpi_next_1=basepair->next_residue(bpi_next);
        string bpj_previous_1=basepair->previous_residue(bpj_previous);
        if(bpj_previous_1!=ResidueConst::RESIDUE_NULL && BPs.count(bpi_next) && basepair->is_basepaired(bpi_next,bpj_previous_1)) {
          hx_c1.insert_residue(bpi_next,Chain::end);
          hx_c2.insert_residue(bpj_previous,Chain::begin);
          hx_c2.insert_residue(bpj_previous_1,Chain::begin);
          bpj_previous=bpj_previous_1;
          BPs.erase(bpi_next);
        }
        else if(bpi_next_1!=ResidueConst::RESIDUE_NULL && BPs.count(bpi_next_1) && basepair->is_basepaired(bpi_next_1,bpj_previous)) {
          hx_c1.insert_residue(bpi_next,Chain::end);
          hx_c1.insert_residue(bpi_next_1,Chain::end);
          hx_c2.insert_residue(bpj_previous,Chain::begin);
          bpi_next=bpi_next_1;
          BPs.erase(bpi_next);
        }
        else if(bpi_next_1!=ResidueConst::RESIDUE_NULL && BPs.count(bpi_next_1) && bpj_previous_1!=ResidueConst::RESIDUE_NULL && basepair->is_basepaired(bpi_next_1,bpj_previous_1)) {
          hx_c1.insert_residue(bpi_next,Chain::end);
          hx_c1.insert_residue(bpi_next_1,Chain::end);
          hx_c2.insert_residue(bpj_previous,Chain::begin);
          hx_c2.insert_residue(bpj_previous_1,Chain::begin);
          bpi_next=bpi_next_1;
          bpj_previous=bpj_previous_1;
          BPs.erase(bpi_next);
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
      //cout<<helices.size()<<":"<<endl;
      //for(auto res : hx.Vresnames)
        //cout<<res<<endl;
      string bp_s=hx_c1.first_of_residue();
      string bp_e=hx_c2.last_of_residue();
      hx_ends.insert(pair<string,string>(bp_s,bp_e));
      bp_s=hx_c2.first_of_residue();
      bp_e=hx_c1.last_of_residue();
      hx_ends.insert(pair<string,string>(bp_s,bp_e));
    }
  }
}

void ExtractRNATemplate::get_loops() {
  loops.clear();
  for(auto& it : hx_ends) {
    string res_i=it.second;
    Chain loop_c(cbs);
    loop_c.insert_residue(res_i,Chain::begin);
    while(1) {
      res_i=basepair->next_residue(res_i);
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
  for(size_t ii=0;ii<basepair->Vresnames.size();ii++) {
    string res_i=basepair->Vresnames[ii];
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
      res_i=basepair->next_residue(res_i);
      if(res_i==ResidueConst::RESIDUE_NULL) break;
      loop_c.insert_residue(res_i,Chain::end);
      if(hx_ends.count(res_i)) break;
      for(size_t jj=0;jj<unsigned_res.size();jj++) {
        if(res_i==unsigned_res[jj]) {
          unsigned_res.erase(unsigned_res.begin()+jj);
          break;
        }
      }
      if(res_i==basepair->Vresnames.back()) break;
    }
    Molecule loop(cbs);
    loop_c.chain_label=Chain::rna;
    loop.insert_chain(loop_c);
    loops.push_back(loop);
  }
}

void ExtractRNATemplate::extract_helix() {
  vector<Molecule> Thelix=helices;
  helices.clear();
  for(size_t ii=0;ii<Thelix.size();ii++) {
    string bp_s=Thelix[ii].chains[0].first_of_residue();
    string bp_e=Thelix[ii].chains[1].last_of_residue();
    //cout<<"extract_helix: "<<bp_s<<" "<<bp_e<<endl;
    int bp_level=basepair->get_basepair_level(bp_s,bp_e);
    if(bp_level==1) {
      RNAMotif hx(cbs);
      hx.type="Helix";
      hx.insert_RNA_segment(Thelix[ii].chains[0],0);
      string secstr_one=basepair->get_segment_secondary_structure(bp_s,Thelix[ii].chains[0].last_of_residue());
      hx.secstrs.push_back(secstr_one);
      hx.insert_RNA_segment(Thelix[ii].chains[1],0);
      secstr_one=basepair->get_segment_secondary_structure(Thelix[ii].chains[1].first_of_residue(),bp_e);
      hx.secstrs.push_back(secstr_one);
      templates.push_back(hx);
    }
    else
      helices.push_back(Thelix[ii]);
  }
}

void ExtractRNATemplate::extract_tail() {
  vector<Molecule> Tloop=loops;
  loops.clear();
  for(size_t ii=0;ii<Tloop.size();ii++) {
    string res_s=Tloop[ii].chains[0].first_of_residue();
    string res_e=Tloop[ii].chains[0].last_of_residue();
    string res_s_p=basepair->previous_residue(res_s);
    string res_e_n=basepair->next_residue(res_e);
    string secstr_one=basepair->get_segment_secondary_structure(res_s,res_e);
    RNAMotif mtf(cbs);
    if((res_s_p==ResidueConst::RESIDUE_NULL && res_e_n==ResidueConst::RESIDUE_NULL) || (res_s_p==res_e && res_e_n==res_s)) {
      mtf.type="Unpaired_chain";
      mtf.insert_RNA_segment(Tloop[ii].chains[0],0);
      mtf.secstrs.push_back(secstr_one);
      //templates.push_back(mtf);
    }
    else if(res_s_p==ResidueConst::RESIDUE_NULL && res_e_n!=ResidueConst::RESIDUE_NULL) {
      mtf.type="Tail5";
      mtf.insert_RNA_segment(Tloop[ii].chains[0],-1);
      mtf.secstrs.push_back(secstr_one);
      templates.push_back(mtf);
    }
    else if(res_s_p!=ResidueConst::RESIDUE_NULL && res_e_n==ResidueConst::RESIDUE_NULL) {
      mtf.type="Tail3";
      mtf.insert_RNA_segment(Tloop[ii].chains[0],-1);
      mtf.secstrs.push_back(secstr_one);
      templates.push_back(mtf);
    }
    else
      loops.push_back(Tloop[ii]);
  }
}

void ExtractRNATemplate::extract_hairpin() {
  vector<Molecule> Tloop=loops;
  loops.clear();
  for(size_t ii=0;ii<Tloop.size();ii++) {
    string res_s=Tloop[ii].chains[0].first_of_residue();
    string res_e=Tloop[ii].chains[0].last_of_residue();
    string secstr_one=basepair->get_segment_secondary_structure(res_s,res_e);
    if(basepair->is_basepaired(res_s,res_e)) {
      RNAMotif hp(cbs);
      hp.type="Hairpin";
      hp.insert_RNA_segment(Tloop[ii].chains[0],-2);
      hp.secstrs.push_back(secstr_one);
      templates.push_back(hp);
    }
    else
      loops.push_back(Tloop[ii]);
  }
}

void ExtractRNATemplate::extract_loop_template() {
  for(size_t ii=0;ii<loops.size();ii++) {
    RNAMotif lp(cbs);
    lp.type="Loop";
    Chain loop_l(cbs);
    loop_l.set_chain(loops[ii].chains[0]);
    lp.insert_RNA_segment(loop_l,-2);
    string res_s_l=loop_l.first_of_residue();
    string res_e_l=loop_l.last_of_residue();
    string secstr_one=basepair->get_segment_secondary_structure(res_s_l,res_e_l);
    lp.secstrs.push_back(secstr_one);
    templates.push_back(lp);
  }
}

void ExtractRNATemplate::extract_pseudoknot_template() {
  map<string,int> hx_end_loop;
  for(auto& it : hx_ends) {
    int nloop=-1;
    for(size_t ii=0;ii<loops.size();ii++) {
      string res_s=loops[ii].chains[0].first_of_residue();
      if(res_s==it.second) {
        nloop=ii;
        break;
      }
    }
    if(nloop>=0) hx_end_loop.insert(pair<string,int>(it.second,nloop));

    nloop=-1;
    for(size_t ii=0;ii<loops.size();ii++) {
      string res_e=loops[ii].chains[0].last_of_residue();
      if(res_e==it.first) {
        nloop=ii;
        break;
      }
    }
    if(nloop>=0) hx_end_loop.insert(pair<string,int>(it.first,nloop));
  }
  vector<int> Vloop;
  vector<int> Vhelix;
  for(size_t ii=0;ii<helices.size();ii++) {
    string S1_a=helices[ii].chains[0].first_of_residue();
    string S1_b=helices[ii].chains[0].last_of_residue();
    string S1_c=helices[ii].chains[1].first_of_residue();
    string S1_d=helices[ii].chains[1].last_of_residue();
    //cout<<"extract_pseudoknot_template: "<<S1_a<<" "<<S1_d<<endl;
    int bp_level=basepair->get_basepair_level(S1_a,S1_d);
    if(bp_level==1) continue;
    if(hx_end_loop.count(S1_b)==0 || hx_end_loop.count(S1_c)==0 || (hx_end_loop.count(S1_a)==0 && hx_end_loop.count(S1_d)==0)) continue;

    for(size_t jj=0;jj<helices.size();jj++) {
      if(ii==jj) continue;
      string S2_a=helices[jj].chains[0].first_of_residue();
      string S2_b=helices[jj].chains[0].last_of_residue();
      string S2_c=helices[jj].chains[1].first_of_residue();
      string S2_d=helices[jj].chains[1].last_of_residue();
      if(hx_end_loop.count(S2_b)==0 || hx_end_loop.count(S2_c)==0 || (hx_end_loop.count(S2_a)==0 && hx_end_loop.count(S2_d)==0)) continue;

      if(hx_end_loop.count(S1_d) && hx_end_loop.count(S2_a)) {
        int L1=hx_end_loop[S1_b];
        int L2=hx_end_loop[S1_c];
        int L3=hx_end_loop[S1_d];
        if(L1==hx_end_loop[S2_a] && L2==hx_end_loop[S2_b] && L3==hx_end_loop[S2_c]) {
          Vloop.push_back(L1);
          Vloop.push_back(L2);
          Vloop.push_back(L3);
          Vhelix.push_back(ii);
          Vhelix.push_back(jj);
          Chain loop_s(cbs), loop_e(cbs);
          loop_s.set_chain(loops[L1].chains[0]);
          loop_e.set_chain(loops[L2].chains[0]);
          RNAMotif mtf(cbs);
          mtf.type="Hairpin_kiss";
          Chain loop_f(cbs);
          loop_f.set_chain(loop_s);
          for(size_t ll=1;ll<helices[jj].chains[0].size()-1;ll++)
            loop_f.insert_residue(helices[jj].chains[0].residues[ll],Chain::end);
          for(size_t ll=0;ll<loop_e.residues.size();ll++)
            loop_f.insert_residue(loop_e.residues[ll],Chain::end);
          mtf.insert_RNA_segment(loop_f,-2);
          mtf.lengths.pop_back();
          int len=(int)loop_s.size()-2;
          mtf.lengths.push_back(len);
          len=(int)helices[jj].chains[0].size();
          mtf.lengths.push_back(len);
          len=(int)loop_e.size()-2;
          mtf.lengths.push_back(len);
          string res_s_f=loop_f.first_of_residue();
          string res_e_f=loop_f.last_of_residue();
          string secstr_one=basepair->get_segment_secondary_structure(res_s_f,res_e_f);
          mtf.secstrs.push_back(secstr_one);
          templates.push_back(mtf);
          loop_s.set_chain(loops[L2].chains[0]);
          loop_e.set_chain(loops[L3].chains[0]);
          mtf.clear();
          mtf.type="Hairpin_kiss";
          loop_f.set_chain(loop_s);
          for(size_t ll=1;ll<helices[ii].chains[1].size()-1;ll++)
            loop_f.insert_residue(helices[ii].chains[1].residues[ll],Chain::end);
          for(size_t ll=0;ll<loop_e.residues.size();ll++)
            loop_f.insert_residue(loop_e.residues[ll],Chain::end);
          mtf.insert_RNA_segment(loop_f,-2);
          mtf.lengths.pop_back();
          len=(int)loop_s.size()-2;
          mtf.lengths.push_back(len);
          len=(int)helices[ii].chains[1].size();
          mtf.lengths.push_back(len);
          len=(int)loop_e.size()-2;
          mtf.lengths.push_back(len);
          res_s_f=loop_f.first_of_residue();
          res_e_f=loop_f.last_of_residue();
          secstr_one=basepair->get_segment_secondary_structure(res_s_f,res_e_f);
          mtf.secstrs.push_back(secstr_one);
          templates.push_back(mtf);
        }
      }

      if(hx_end_loop.count(S1_a) && hx_end_loop.count(S2_d)) {
        int L1=hx_end_loop[S2_b];
        int L2=hx_end_loop[S2_c];
        int L3=hx_end_loop[S2_d];
        if(L1==hx_end_loop[S1_a] && L2==hx_end_loop[S1_b] && L3==hx_end_loop[S1_c]) {
          Vloop.push_back(L1);
          Vloop.push_back(L2);
          Vloop.push_back(L3);
          Vhelix.push_back(ii);
          Vhelix.push_back(jj);
          Chain loop_s(cbs), loop_e(cbs);
          loop_s.set_chain(loops[L1].chains[0]);
          loop_e.set_chain(loops[L2].chains[0]);
          RNAMotif mtf(cbs);
          mtf.type="Hairpin_kiss";
          Chain loop_f(cbs);
          loop_f.set_chain(loop_s);
          for(size_t ll=1;ll<helices[ii].chains[0].size()-1;ll++)
            loop_f.insert_residue(helices[ii].chains[0].residues[ll],Chain::end);
          for(size_t ll=0;ll<loop_e.residues.size();ll++)
            loop_f.insert_residue(loop_e.residues[ll],Chain::end);
          mtf.insert_RNA_segment(loop_f,-2);
          mtf.lengths.pop_back();
          int len=(int)loop_s.size()-2;
          mtf.lengths.push_back(len);
          len=(int)helices[ii].chains[0].size();
          mtf.lengths.push_back(len);
          len=(int)loop_e.size()-2;
          mtf.lengths.push_back(len);
          string res_s_f=loop_f.first_of_residue();
          string res_e_f=loop_f.last_of_residue();
          string secstr_one=basepair->get_segment_secondary_structure(res_s_f,res_e_f);
          mtf.secstrs.push_back(secstr_one);
          templates.push_back(mtf);
          loop_s.set_chain(loops[L2].chains[0]);
          loop_e.set_chain(loops[L3].chains[0]);
          mtf.clear();
          mtf.type="Hairpin_kiss";
          loop_f.set_chain(loop_s);
          for(size_t ll=1;ll<helices[jj].chains[1].size()-1;ll++)
            loop_f.insert_residue(helices[jj].chains[1].residues[ll],Chain::end);
          for(size_t ll=0;ll<loop_e.residues.size();ll++)
            loop_f.insert_residue(loop_e.residues[ll],Chain::end);
          mtf.insert_RNA_segment(loop_f,-2);
          mtf.lengths.pop_back();
          len=(int)loop_s.size()-2;
          mtf.lengths.push_back(len);
          len=(int)helices[jj].chains[1].size();
          mtf.lengths.push_back(len);
          len=(int)loop_e.size()-2;
          mtf.lengths.push_back(len);
          res_s_f=loop_f.first_of_residue();
          res_e_f=loop_f.last_of_residue();
          secstr_one=basepair->get_segment_secondary_structure(res_s_f,res_e_f);
          mtf.secstrs.push_back(secstr_one);
          templates.push_back(mtf);
        }
      }
    }
  }

  vector<Molecule> Thelix=helices;
  helices.clear();
  for(size_t ii=0;ii<Thelix.size();ii++) {
    if(count(Vhelix.begin(),Vhelix.end(),ii)) {
      string bp_s=Thelix[ii].chains[0].first_of_residue();
      hx_ends.erase(bp_s);
      bp_s=Thelix[ii].chains[1].first_of_residue();
      hx_ends.erase(bp_s);
    }
    else helices.push_back(Thelix[ii]);
  }

  vector<Molecule> Tloop=loops;
  loops.clear();
  for(size_t ii=0;ii<Tloop.size();ii++) {
    if(count(Vloop.begin(),Vloop.end(),ii)) continue;
    else loops.push_back(Tloop[ii]);
  }
}

void ExtractRNATemplate::extract_kiss_template() {
  vector<int> Vloop;
  for(size_t ii=0;ii<helices.size();ii++) {
    for(int jj=0;jj<2;jj++) {
      string res_s=helices[ii].chains[jj].first_of_residue();
      string res_e=helices[ii].chains[jj].last_of_residue();
      //cout<<"extract_kiss_template: "<<res_s<<" "<<hx_ends[res_s]<<endl;
      int bp_level=basepair->get_basepair_level(res_s,hx_ends[res_s]);
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
          Vloop.push_back(kk);
        }
      }
      for(size_t kk=0;kk<loops.size();kk++) {
        string res_s_l=loops[kk].chains[0].first_of_residue();
        if(res_s_l==res_e) {
          check_e=1;
          loop_e.set_chain(loops[kk].chains[0]);
          Vloop.push_back(kk);
        }
      }

      if(check_s && check_e) {
        RNAMotif mtf(cbs);
        string res_s_m=loop_s.first_of_residue();
        string res_e_m=loop_e.last_of_residue();
        if(basepair->is_basepaired(res_s_m,res_e_m)) mtf.type="Hairpin_kiss";
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
        string secstr_one=basepair->get_segment_secondary_structure(res_s_f,res_e_f);
        mtf.secstrs.push_back(secstr_one);
        templates.push_back(mtf);
      }
    }
  }

  vector<Molecule> Tloop=loops;
  loops.clear();
  for(size_t ii=0;ii<Tloop.size();ii++) {
    if(count(Vloop.begin(),Vloop.end(),ii)) continue;
    else loops.push_back(Tloop[ii]);
  }
}

void ExtractRNATemplate::extract_junction_template() {
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
      int ns=0;
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
        Chain loop_l(cbs);
        loop_l.set_chain(loops[Vloop[jj]].chains[0]);
        mtf.insert_RNA_segment(loop_l,-2);
        string res_s_l=loop_l.first_of_residue();
        string res_e_l=loop_l.last_of_residue();
        string secstr_one=basepair->get_segment_secondary_structure(res_s_l,res_e_l);
        mtf.secstrs.push_back(secstr_one);
      }
      for(int jj=0;jj<ns;jj++) {
        Chain loop_l(cbs);
        loop_l.set_chain(loops[Vloop[jj]].chains[0]);
        mtf.insert_RNA_segment(loop_l,-2);
        string res_s_l=loop_l.first_of_residue();
        string res_e_l=loop_l.last_of_residue();
        string secstr_one=basepair->get_segment_secondary_structure(res_s_l,res_e_l);
        mtf.secstrs.push_back(secstr_one);
      }

      vector<RNAMotif> Vmotif=rotate_junction(mtf);
      for(size_t kk=0;kk<Vmotif.size();kk++)
        templates.push_back(Vmotif[kk]);

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
      res_s_m=basepair->previous_residue(res_s_m);
      string res_e_m=loops[Vloop.back()].chains[0].last_of_residue();
      res_e_m=basepair->next_residue(res_e_m);
      for(size_t jj=0;jj<templates.size();jj++) {
        if(templates[jj].type!="Loop_kiss_loop") continue;
        string res_s_k=templates[jj].chains[0].first_of_residue();
        string res_e_k=templates[jj].chains[0].last_of_residue();
        if(basepair->is_basepaired(res_s_k,res_e_m) && basepair->is_basepaired(res_e_k,res_s_m)) {
          RNAMotif mtf(cbs);
          mtf.type=to_string((int)Vloop.size()+1)+"-way_junction_kiss";
          for(size_t kk=0;kk<Vloop.size();kk++) {
            Chain loop_l(cbs);
            loop_l.set_chain(loops[Vloop[kk]].chains[0]);
            mtf.insert_RNA_segment(loop_l,-2);
            string res_s_l=loop_l.first_of_residue();
            string res_e_l=loop_l.last_of_residue();
            string secstr_one=basepair->get_segment_secondary_structure(res_s_l,res_e_l);
            mtf.secstrs.push_back(secstr_one);
          }
          mtf.insert_RNA_segment(templates[jj].chains[0],-2);
          mtf.secstrs.push_back(templates[jj].secstrs[0]);
          templates.push_back(mtf);

          vector<Molecule> Tloop=loops;
          loops.clear();
          for(size_t jj=0;jj<Tloop.size();jj++) {
            if(count(Vloop.begin(),Vloop.end(),jj)) continue;
            else loops.push_back(Tloop[jj]);
          }
        }
      }
    }
  }
}

vector<RNAMotif> ExtractRNATemplate::rotate_junction(RNAMotif &motif_one) {
  vector<RNAMotif> Vmotif;
  int len_min=motif_one.lengths[0];
  for(size_t ii=1;ii<motif_one.size();ii++)
    if(motif_one.lengths[ii]<=len_min) 
      len_min=motif_one.lengths[ii];

  for(size_t ii=0;ii<motif_one.size();ii++) {
    if(motif_one.lengths[ii]==len_min) {
      RNAMotif mtf(cbs);
      for(size_t jj=ii;jj<motif_one.size();jj++) {
        mtf.insert_RNA_segment(motif_one.chains[jj],-2);
        string secstr_one=motif_one.secstrs[jj];
        mtf.secstrs.push_back(secstr_one);
      }
      for(size_t jj=0;jj<ii;jj++) {
        mtf.insert_RNA_segment(motif_one.chains[jj],-2);
        string secstr_one=motif_one.secstrs[jj];
        mtf.secstrs.push_back(secstr_one);
      }
      mtf.type=motif_one.type;
      Vmotif.push_back(mtf);
    }
  }

  return Vmotif;
}

inline bool ExtractRNATemplate::is_pseudo_basepairs(string &nt_i) {
  if(pseudo_basepairs.count(nt_i)==0) return 0;
  else return 1;
}

inline bool ExtractRNATemplate::is_pseudo_basepairs(string &nt_i, string &nt_j) {
  if(pseudo_basepairs.count(nt_i)==0) return 0;
  else {
    if(nt_j==pseudo_basepairs[nt_i]) return 1;
    else return 0;
  }
}

void ExtractRNATemplate::extend_loop_kiss_loop(RNAMotif &motif_one) {
  vector<int> Vend1,Vend2,Vend3,Vend4;
  string nt_previous=basepair->previous_residue(motif_one.chains[0].first_of_residue());
  if(nt_previous!=ResidueConst::RESIDUE_NULL) Vend1.push_back(-1);
  Vend1.push_back(0);
  string nt_one=motif_one.chains[0][0].get_residue_name();
  if(ExtractRNATemplate::is_pseudo_basepairs(nt_one)) Vend1.push_back(1);
  int len1=motif_one.lengths[0];
  int len2=motif_one.lengths[1];
  int len3=motif_one.lengths[2];
  nt_one=motif_one.chains[0][len1].get_residue_name();
  if(ExtractRNATemplate::is_pseudo_basepairs(nt_one)) Vend2.push_back(-1);
  Vend2.push_back(0);
  Vend2.push_back(1);
  Vend3.push_back(-1);
  Vend3.push_back(0);
  nt_one=motif_one.chains[0][len1+len2+1].get_residue_name();
  if(ExtractRNATemplate::is_pseudo_basepairs(nt_one)) Vend3.push_back(1);
  nt_one=motif_one.chains[0][len1+len2+len3].get_residue_name();
  if(ExtractRNATemplate::is_pseudo_basepairs(nt_one)) Vend4.push_back(-1);
  Vend4.push_back(0);
  string nt_next=basepair->next_residue(motif_one.chains[0].last_of_residue());
  if(nt_next!=ResidueConst::RESIDUE_NULL) Vend4.push_back(1);
  
  for(size_t e1=0;e1<Vend1.size();e1++) {
    for(size_t e2=0;e2<Vend2.size();e2++) {
      for(size_t e3=0;e3<Vend3.size();e3++) {
        for(size_t e4=0;e4<Vend4.size();e4++) {
          if(Vend1[e1]==0 && Vend2[e2]==0 && Vend3[e3]==0 && Vend4[e4]==0) continue;
          Chain LKL=motif_one.chains[0];
          if(Vend1[e1]==-1)
            LKL.insert_residue(nt_previous,Chain::begin);
          else if(Vend1[e1]==1)
            LKL.delete_residue(Chain::begin);
          if(Vend4[e4]==-1)
            LKL.delete_residue(Chain::end);
          else if(Vend4[e4]==1)
            LKL.insert_residue(nt_next,Chain::end);
          int l1=len1-Vend1[e1]+Vend2[e2];
          int l2=len2-Vend2[e2]+Vend3[e3];
          int l3=len3-Vend3[e3]+Vend4[e4];
          if(l1<0 || l3<0 || l2<2) continue;
          RNAMotif mtf(cbs);
          mtf.type=motif_one.type;
          mtf.insert_RNA_segment(LKL,-2);
          mtf.lengths.pop_back();
          mtf.lengths.push_back(l1);
          mtf.lengths.push_back(l2);
          mtf.lengths.push_back(l3);
          string res_s=LKL.first_of_residue();
          string res_e=LKL.last_of_residue();
          string secstr_one=basepair->get_segment_secondary_structure(res_s,res_e);
          mtf.secstrs.push_back(secstr_one);
          extended_templates.push_back(mtf);
        }
      }
    }
  }

}

void ExtractRNATemplate::extend_hairpin_kiss(RNAMotif &motif_one) {
  int len1=motif_one.lengths[0];
  int len2=motif_one.lengths[1];
  int len3=motif_one.lengths[2];
  vector<int> Vend1,Vend2,Vend3;
  string nt_previous=basepair->previous_residue(motif_one.chains[0].first_of_residue());
  string nt_next=basepair->next_residue(motif_one.chains[0].last_of_residue());
  if(nt_previous!=ResidueConst::RESIDUE_NULL && nt_next!=ResidueConst::RESIDUE_NULL) 
    Vend1.push_back(-1);
  Vend1.push_back(0);
  string nt_i=motif_one.chains[0][1].get_residue_name();
  string nt_j=motif_one.chains[0][len1+len2+len3].get_residue_name();
  if(ExtractRNATemplate::is_pseudo_basepairs(nt_i,nt_j)) Vend1.push_back(1);
  string nt_one=motif_one.chains[0][len1].get_residue_name();
  if(ExtractRNATemplate::is_pseudo_basepairs(nt_one)) Vend2.push_back(-1);
  Vend2.push_back(0);
  Vend2.push_back(1);
  Vend3.push_back(-1);
  Vend3.push_back(0);
  nt_one=motif_one.chains[0][len1+len2+1].get_residue_name();
  if(ExtractRNATemplate::is_pseudo_basepairs(nt_one)) Vend3.push_back(1);
  
  for(size_t e1=0;e1<Vend1.size();e1++) {
    for(size_t e2=0;e2<Vend2.size();e2++) {
      for(size_t e3=0;e3<Vend3.size();e3++) {
        if(Vend1[e1]==0 && Vend2[e2]==0 && Vend3[e3]==0) continue;
        Chain HK=motif_one.chains[0];
        if(Vend1[e1]==-1) {
          HK.insert_residue(nt_previous,Chain::begin);
          HK.insert_residue(nt_next,Chain::end);
        }
        else if(Vend1[e1]==1) {
          HK.delete_residue(Chain::begin);
          HK.delete_residue(Chain::end);
        }
        int l1=len1-Vend1[e1]+Vend2[e2];
        int l2=len2-Vend2[e2]+Vend3[e3];
        int l3=len3-Vend3[e3]-Vend1[e1];
        if(l1<0 || l3<0 || l2<2) continue;
        RNAMotif mtf(cbs);
        mtf.type=motif_one.type;
        mtf.insert_RNA_segment(HK,-2);
        mtf.lengths.pop_back();
        mtf.lengths.push_back(l1);
        mtf.lengths.push_back(l2);
        mtf.lengths.push_back(l3);
        string res_s=HK.first_of_residue();
        string res_e=HK.last_of_residue();
        string secstr_one=basepair->get_segment_secondary_structure(res_s,res_e);
        mtf.secstrs.push_back(secstr_one);
        extended_templates.push_back(mtf);
      }
    }
  }
}

void ExtractRNATemplate::extend_junction(RNAMotif &motif_one) {
  if(motif_one.type=="2-way_junction" && motif_one.lengths[0]<=7 && motif_one.lengths[1]<=7) 
    return;
  vector<string> stoken=split_string(motif_one.type,"-way_");
  int Nw=stoi(stoken[0]);
  if(Nw>=Nw_extend_max) return;
  vector<vector<int> > Vends;
  vector<string> Vnt_previous,Vnt_next;
  string nt_next_last;
  for(int ii=0;ii<Nw;ii++) {
    vector<int> Vend;
    Chain loop_i(cbs),loop_j(cbs);
    loop_i.set_chain(motif_one.chains[ii]);
    if(ii==0) loop_j.set_chain(motif_one.chains.back());
    else loop_j.set_chain(motif_one.chains[ii-1]);
    string nt_previous=basepair->previous_residue(loop_i.first_of_residue());
    string nt_next=basepair->next_residue(loop_j.last_of_residue());
    Vnt_previous.push_back(nt_previous);
    if(ii==0) nt_next_last=nt_next;
    else Vnt_next.push_back(nt_next);
    if(nt_previous!=ResidueConst::RESIDUE_NULL && nt_next!=ResidueConst::RESIDUE_NULL)
      Vend.push_back(-1);
    Vend.push_back(0);
    string nt_i=loop_i[1].get_residue_name();
    string nt_j=loop_j[loop_j.size()-2].get_residue_name();
    if(ExtractRNATemplate::is_pseudo_basepairs(nt_i,nt_j)) Vend.push_back(1);
    Vends.push_back(Vend);
  }
  Vnt_next.push_back(nt_next_last);

  vector<vector<int> > OPs;
  vector<int> Vacc;
  for(int ii=0;ii<Nw;ii++) {
    int n0=1;
    for(int jj=ii+1;jj<Nw;jj++)
      n0*=Vends[jj].size();
    Vacc.push_back(n0);
  }
  int n=0;
  while(n<Vacc[0]*(int)Vends[0].size()) {
    vector<int> op;
    for(int jj=0;jj<Nw;jj++) {
      int index=((int)(n/Vacc[jj]))%Vends[jj].size();
      op.push_back(Vends[jj][index]);
    }
    n++;
    OPs.push_back(op);
  }

  vector<RNAMotif> Vmotif;
  for(size_t ii=0;ii<OPs.size();ii++) {
    RNAMotif mtf(cbs);
    mtf.type=motif_one.type;
    int check1=0;
    int check2=0;
    for(int jj=0;jj<Nw;jj++) {
      Chain loop_i=motif_one.chains[jj];
      int nnext=jj+1;
      if(nnext==Nw) nnext=0;
      if(OPs[ii][jj]==-1)
        loop_i.insert_residue(Vnt_previous[jj],Chain::begin);
      else if(OPs[ii][jj]==0)
        check2++;
      else if(OPs[ii][jj]==1)
        loop_i.delete_residue(Chain::begin);
      if(OPs[ii][nnext]==-1)
        loop_i.insert_residue(Vnt_next[jj],Chain::end);
      else if(OPs[ii][nnext]==1)
        loop_i.delete_residue(Chain::end);
      if(loop_i.size()<2) {
        check1=1; break;
      }
      mtf.insert_RNA_segment(loop_i,-2);
      string res_s=loop_i.first_of_residue();
      string res_e=loop_i.last_of_residue();
      string secstr_one=basepair->get_segment_secondary_structure(res_s,res_e);
      mtf.secstrs.push_back(secstr_one);
      loop_i.clear();
    }
    if(check1) continue;
    if(check2==Nw) continue;

    Vmotif.push_back(mtf);
    mtf.clear();
  }

  if(motif_one.type.find("_kiss")==string::npos) {
    for(size_t jj=0;jj<Vmotif.size();jj++) {
      if(motif_one.type=="2-way_junction" && Vmotif[jj].lengths[0]==0 && Vmotif[jj].lengths[1]==0) continue;
      vector<RNAMotif> Tmotif=rotate_junction(Vmotif[jj]);
      extended_templates.insert(extended_templates.end(),Tmotif.begin(),Tmotif.end());
    }
  }
  else {
    //treat N-way_junction_kiss as N-way_junction in extended templates
    Vmotif.insert(Vmotif.begin(),motif_one);
    for(size_t jj=0;jj<Vmotif.size();jj++) {
      Vmotif[jj].type=to_string(Nw)+"-way_junction";
      vector<RNAMotif> Tmotif=rotate_junction(Vmotif[jj]);
      extended_templates.insert(extended_templates.end(),Tmotif.begin(),Tmotif.end());
    }
  }
}

ostream& operator <<(ostream &output, ExtractRNATemplate &ERtemplate_one) {
  output<<"List of "<<ERtemplate_one.size()<<" templates:"<<endl;
  for(size_t ii=0;ii<ERtemplate_one.size();ii++) {
    string mtf_name=ERtemplate_one.templates[ii].get_RNA_motif_info();
    output<<mtf_name<<endl;
  }
  size_t num=ERtemplate_one.extended_templates.size();
  if(num>0) {
    output<<"---------------------------------------"<<endl;
    output<<"List of "<<num<<" extended templates:"<<endl;
    for(size_t ii=0;ii<num;ii++) {
      string mtf_name=ERtemplate_one.extended_templates[ii].get_RNA_motif_info();
      output<<mtf_name<<endl;
    }
  }

  return output;
}