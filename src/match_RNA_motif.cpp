#include "basic.h"
#include "initialization.h"
#include "error.h"
#include "match_RNA_motif.h"

using namespace CBS_NS;

const int MatchRNAMotif::MATCH_NULL = 100;

MatchRNAMotif::MatchRNAMotif(Cobmms *cbs) : Pointers(cbs) {
  match_score = MATCH_NULL;
  motif_ptr = nullptr;
  template_ptr = nullptr;  
}

MatchRNAMotif::MatchRNAMotif(const MatchRNAMotif &other) : Pointers(other.cbs) {
  match_score=other.match_score;
  //length_difference=other.length_difference;
  sequence_similarity=other.sequence_similarity;
  match_code=other.match_code;
  template_info=other.template_info;

  motif_ptr = nullptr;
  template_ptr = nullptr;
}

MatchRNAMotif& MatchRNAMotif::operator=(const MatchRNAMotif &other) {
  if(this != &other) {
    Pointers::operator=(other.cbs);
    match_score=other.match_score;
    //length_difference=other.length_difference;
    sequence_similarity=other.sequence_similarity;
    match_code=other.match_code;
    template_info=other.template_info;
  
    motif_ptr = nullptr;
    template_ptr = nullptr;
  }

  return *this;
}

MatchRNAMotif::~MatchRNAMotif() {
  match_score=MATCH_NULL;
  //length_difference=MATCH_NULL;
  sequence_similarity=0.0;
  match_code="";
  template_info="";  

  motif_ptr = nullptr;
  template_ptr = nullptr; 
}


void MatchRNAMotif::get_match_score(const RNAMotif &motif_one, string &template_one) {
  motif_ptr=&motif_one;
  template_info=template_one;
  template_ptr= new RNAMotif(cbs);
  template_ptr->set_RNA_motif(template_one,1);

  match_code="";

  if(motif_ptr->type!=template_ptr->type) {
    string str="match_RNA_motif.cpp: Type of motif is not matched --> "+motif_ptr->type+"/"+template_ptr->type;
    error->one(str);  
  }
  if(motif_ptr->lengths.size()!=template_ptr->lengths.size()) {
    string str="match_RNA_motif.cpp: Name of motif is not matched --> ";
    str+=motif_ptr->get_RNA_motif_name()+"/"+template_ptr->get_RNA_motif_name();
    error->one(str);
  }

  if(motif_ptr->type=="Loop") MatchRNAMotif::match_Loop();
  else if(motif_ptr->type=="Hairpin") MatchRNAMotif::match_Hairpin();
  else if(motif_ptr->type=="Tail5" || motif_ptr->type=="Tail3") MatchRNAMotif::match_Tail();
  else if(motif_ptr->type=="Hairpin_kiss" || motif_ptr->type=="Loop_kiss_loop") MatchRNAMotif::match_Loop_kiss();
  else if(motif_ptr->type.find("way_junction")!=string::npos) MatchRNAMotif::match_Junction();
  else {
    string str="match_RNA_motif.cpp: Unknown motif type --> "+motif_ptr->type;
    error->one(str);  
  }

  delete template_ptr;
  template_ptr = nullptr;
  motif_ptr = nullptr;
}

string MatchRNAMotif::get_match_info() const {
  string str=match_code+" "+to_string(sequence_similarity)+" "+template_info;
  
  return str;
}

void MatchRNAMotif::clear() {
  match_score=MATCH_NULL;
  //length_difference=MATCH_NULL;
  sequence_similarity=0.0;
  match_code="";
  template_info=""; 

  motif_ptr = nullptr;
  template_ptr = nullptr; 
}

void MatchRNAMotif::match_Loop() {
  if(motif_ptr->lengths[0]==template_ptr->lengths[0]) {
    match_score=0;
    match_code="1";
  }  
  else {
    match_score=MATCH_NULL;
    match_code="0";
  }
  //length_difference=abs(motif_ptr->lengths[0]-template_ptr->lengths[0])*1.0/(motif_ptr->lengths[0]+2.0);
  string seq_ref=motif_ptr->sequences[0];
  int flag_check=1;
  if(motif_ptr->lengths[0] <= 5)
    flag_check=0;
  //seq_ref=seq_ref.substr(1,seq_ref.size()-2);
  if(flag_check) {
    MatchRNAMotif::modify_paired_nucleotide(seq_ref,0);
    MatchRNAMotif::modify_paired_nucleotide(seq_ref,seq_ref.size()-1);
  }
  string seq_tar=template_ptr->sequences[0];
  //seq_tar=seq_tar.substr(1,seq_tar.size()-2);
  if(flag_check) {
    MatchRNAMotif::modify_paired_nucleotide(seq_tar,0);
    MatchRNAMotif::modify_paired_nucleotide(seq_tar,seq_tar.size()-1);
  }
  sequence_similarity=MatchRNAMotif::calculate_sequence_similarity(seq_ref,seq_tar);
}

void MatchRNAMotif::match_Hairpin() {
  if(motif_ptr->lengths[0]==template_ptr->lengths[0]) {
    match_score=0;
    match_code="1";
  }  
  else {
    match_score=MATCH_NULL;
    match_code="0";
  }
  //length_difference=abs(motif_ptr->lengths[0]-template_ptr->lengths[0])*1.0/(motif_ptr->lengths[0]+2.0); 
  string seq_ref=motif_ptr->sequences[0];
  int flag_check=1;
  if(motif_ptr->lengths[0] >= 3 && motif_ptr->lengths[0] <= 6)
    flag_check=0;
  //seq_ref=seq_ref.substr(1,seq_ref.size()-2);
  if(flag_check) {
    MatchRNAMotif::modify_paired_nucleotide(seq_ref,0);
    MatchRNAMotif::modify_paired_nucleotide(seq_ref,seq_ref.size()-1);
  }
  string seq_tar=template_ptr->sequences[0];
  //seq_tar=seq_tar.substr(1,seq_tar.size()-2);
  if(flag_check) {
    MatchRNAMotif::modify_paired_nucleotide(seq_tar,0);
    MatchRNAMotif::modify_paired_nucleotide(seq_tar,seq_tar.size()-1);
  }
  sequence_similarity=MatchRNAMotif::calculate_sequence_similarity(seq_ref,seq_tar);
}

void MatchRNAMotif::match_Tail() {
  if(motif_ptr->lengths[0]==template_ptr->lengths[0]) {
    match_score=0;
    match_code="1";
  }
  else {
    match_score=MATCH_NULL;
    match_code="0";
  }
  //length_difference=abs(motif_ptr->lengths[0]-template_ptr->lengths[0])*1.0/(motif_ptr->lengths[0]+1.0); 
  string seq_ref=motif_ptr->sequences[0];
  string seq_tar=template_ptr->sequences[0];
  int flag_check=1;
  if(motif_ptr->lengths[0] <= 5)
    flag_check=0;
  if(motif_ptr->type=="Tail5") {
    //seq_ref=seq_ref.substr(0,seq_ref.size()-1);
    //seq_tar=seq_tar.substr(0,seq_tar.size()-1); 
    if(flag_check) {
      MatchRNAMotif::modify_paired_nucleotide(seq_ref,seq_ref.size()-1);
      MatchRNAMotif::modify_paired_nucleotide(seq_tar,seq_tar.size()-1);
    }
  }
  else {
    //seq_ref=seq_ref.substr(1);
    //seq_tar=seq_tar.substr(1);
    if(flag_check) {
      MatchRNAMotif::modify_paired_nucleotide(seq_ref,0);
      MatchRNAMotif::modify_paired_nucleotide(seq_tar,0);  
    }
  }
  sequence_similarity=MatchRNAMotif::calculate_sequence_similarity(seq_ref,seq_tar);
}

void MatchRNAMotif::match_Loop_kiss() {
  match_score=0;
  for(size_t ii=0;ii<3;ii++) {
    if(motif_ptr->lengths[ii]==template_ptr->lengths[ii]) 
      match_code+="1";
    else {
      match_score+=1;
      match_code+="0";
    }
  }
  if(match_code=="000") match_score=MATCH_NULL;
  else if(match_code=="010") {
    match_code="000";
    match_score=MATCH_NULL;
  }
  else if(match_code=="101") {
    if(motif_ptr->lengths[0] <= motif_ptr->lengths[2])
      match_code="100";
    else
      match_code="001";
    match_score=2;
  }
  //length_difference=abs(motif_ptr->lengths[0]-template_ptr->lengths[0])*1.0/(motif_ptr->lengths[0]+1.0);
  //length_difference+=abs(motif_ptr->lengths[1]-template_ptr->lengths[1])*1.0/motif_ptr->lengths[1];
  //length_difference+=abs(motif_ptr->lengths[2]-template_ptr->lengths[2])*1.0/(motif_ptr->lengths[2]+1.0);
  //string seq_ref=motif_ptr->sequences[0].substr(1,motif_ptr->lengths[0]);
  //int index_ref=motif_ptr->lengths[0]+motif_ptr->lengths[1]+1;
  //seq_ref+=" "+motif_ptr->sequences[0].substr(index_ref,motif_ptr->lengths[2]);
  int flag_check=1;
  if(motif_ptr->lengths[0] <= 1 && motif_ptr->lengths[2] <= 1) {
    flag_check=0;
  }
  string seq_ref=motif_ptr->sequences[0];
  if(flag_check)
    MatchRNAMotif::modify_paired_nucleotide(seq_ref,0);
  for(int ii=0;ii < motif_ptr->lengths[1];ii++) {
    if(ii==0 && flag_check==0) continue;
    else if(ii==motif_ptr->lengths[1]-1 && flag_check==0) continue;
    size_t index_ref=motif_ptr->lengths[0]+ii+1;
    MatchRNAMotif::modify_paired_nucleotide(seq_ref,index_ref);
  }
  if(flag_check)
    MatchRNAMotif::modify_paired_nucleotide(seq_ref,seq_ref.size()-1);
  //string seq_tar=template_ptr->sequences[0].substr(1,template_ptr->lengths[0]);
  //int index_tar=template_ptr->lengths[0]+template_ptr->lengths[1]+1;
  //seq_tar+=" "+template_ptr->sequences[0].substr(index_tar,template_ptr->lengths[2]);
  string seq_tar=template_ptr->sequences[0];
  if(flag_check)
    MatchRNAMotif::modify_paired_nucleotide(seq_tar,0);
  for(int ii=0;ii < template_ptr->lengths[1];ii++) {
    if(ii==0 && flag_check==0) continue;
    else if(ii==template_ptr->lengths[1]-1 && flag_check==0) continue;
    size_t index_tar=template_ptr->lengths[0]+ii+1;
    MatchRNAMotif::modify_paired_nucleotide(seq_tar,index_tar);
  }
  if(flag_check)
    MatchRNAMotif::modify_paired_nucleotide(seq_tar,seq_tar.size()-1);
  sequence_similarity=MatchRNAMotif::calculate_sequence_similarity(seq_ref,seq_tar);
}

void MatchRNAMotif::match_Junction() {
  match_score=0;
  //length_difference=0.0;
  int len_total=0;
  for(size_t ii=0;ii<motif_ptr->lengths.size();ii++) {
    len_total+=motif_ptr->lengths[ii];
  }
  int flag_check=1;
  if(motif_ptr->type=="2-way_junction" && motif_ptr->lengths[0] < 5 && motif_ptr->lengths[1] < 5)
    flag_check=0;
  if(len_total<2)
    flag_check=0;
  string seq_ref,seq_tar;
  for(size_t ii=0;ii<motif_ptr->lengths.size();ii++) {
    if(motif_ptr->lengths[ii]==template_ptr->lengths[ii]) {
      match_code+="1";
    }
    else {
      match_score+=1;
      match_code+="0";
    }
    //length_difference+=abs(motif_ptr->lengths[ii]-template_ptr->lengths[ii])*1.0/(motif_ptr->lengths[ii]+2.0);
    string seq_ref_i=motif_ptr->sequences[ii];
    //seq_ref_i=seq_ref_i.substr(1,seq_ref_i.size()-2);
    if(flag_check) {
      MatchRNAMotif::modify_paired_nucleotide(seq_ref_i,0);
      MatchRNAMotif::modify_paired_nucleotide(seq_ref_i,seq_ref_i.size()-1);
    }
    string seq_tar_i=template_ptr->sequences[ii];
    //seq_tar_i=seq_tar_i.substr(1,seq_tar_i.size()-2);
    if(flag_check) {
      MatchRNAMotif::modify_paired_nucleotide(seq_tar_i,0);
      MatchRNAMotif::modify_paired_nucleotide(seq_tar_i,seq_tar_i.size()-1);
    }
    if(seq_ref_i.size()<seq_tar_i.size()) {
      size_t len_diff=seq_tar_i.size()-seq_ref_i.size();
      size_t index_half=round(seq_ref_i.size()/2);
      string seq_h2=seq_ref_i.substr(index_half);
      seq_ref_i=seq_ref_i.substr(0,index_half);
      for(size_t kk=0;kk<len_diff;kk++)
        seq_ref_i+=" ";
      seq_ref_i+=seq_h2;
    }
    else if(seq_ref_i.size()>seq_tar_i.size()) {
      size_t len_diff=seq_ref_i.size()-seq_tar_i.size();
      size_t index_half=round(seq_tar_i.size()/2);
      string seq_h2=seq_tar_i.substr(index_half);
      seq_tar_i=seq_tar_i.substr(0,index_half);
      for(size_t kk=0;kk<len_diff;kk++)
        seq_tar_i+=" ";
      seq_tar_i+=seq_h2;
    }
    if(ii==0) {
      seq_ref=seq_ref_i;
      seq_tar=seq_tar_i;
    }
    else {
      seq_ref+=" "+seq_ref_i;
      seq_tar+=" "+seq_tar_i;
    }
  }
  sequence_similarity=MatchRNAMotif::calculate_sequence_similarity(seq_ref,seq_tar);
  if(match_score==(int)motif_ptr->lengths.size()) 
    match_score=MATCH_NULL;
}

void MatchRNAMotif::modify_paired_nucleotide(string &seq_one, size_t index_one) {
  if(index_one >= seq_one.size()) {
    string str="match_RNA_motif.cpp: Wrong sequence index --> "+to_string(index_one)+" / "+to_string(seq_one.size());
    error->one(str);
  }
  // P = paired A,G,R; Q = paired C,U,Y
  if(seq_one[index_one]=='A' || seq_one[index_one]=='G' || seq_one[index_one]=='R') seq_one[index_one]='P';
  else if(seq_one[index_one]=='C' || seq_one[index_one]=='U' || seq_one[index_one]=='Y') seq_one[index_one]='Q';
}

double MatchRNAMotif::calculate_sequence_similarity(const string &seq1, const string &seq2) {
  if(seq1.size()!=seq2.size()) {
    string str="match_RNA_motif.cpp: Sequence lengths are not matched --> \n"+seq1+"("+to_string(seq1.size())+")";
    str+="\n"+seq2+"("+to_string(seq2.size())+")";
    error->warning(str);
    return 0.0;
  }

  if(seq1.size()==0) return 1.0;

  double score=0.0;
  for(size_t ii=0;ii<seq1.size();ii++) {
    char c1=seq1[ii];
    char c2=seq2[ii];
    if(c1==c2 || (c1=='N' && c2!=' ') || (c1!=' ' && c2=='N' )) score+=1.0;
    else if((c1=='A' && c2=='G') || (c1=='G' && c2=='A')) score+=0.5;
    else if((c1=='C' && c2=='U') || (c1=='U' && c2=='C')) score+=0.5;
    else if((c1=='A' || c1=='G') && (c2=='U' || c2=='C')) score+=0;
    else if((c1=='C' || c1=='U') && (c2=='A' || c2=='G')) score+=0;
    else if((c1=='A' && c2=='R') || (c1=='R' && c2=='A')) score+=1.0;
    else if((c1=='G' && c2=='R') || (c1=='R' && c2=='G')) score+=1.0;
    else if((c1=='C' && c2=='Y') || (c1=='Y' && c2=='C')) score+=1.0;
    else if((c1=='U' && c2=='Y') || (c1=='Y' && c2=='U')) score+=1.0;
    else if((c1=='R' && c2=='Y') || (c1=='Y' && c2=='R')) score+=0;
    else if(c1==' ' && (c2=='A' || c2=='G' || c2=='C' || c2=='U')) score+=-1.0;
    else if(c1==' ' && (c2=='R' || c2=='Y' || c2=='N')) score+=-1.0;
    else if((c1=='A' || c1=='G' || c1=='C' || c1=='U') && c2==' ') score+=-1.0;
    else if((c1=='R' || c1=='Y' || c1=='N') && c2==' ') score+=-1.0;
    else if((c1=='P' && c2=='Q') || (c1=='Q' && c2=='P')) score+=0.5;
    else if(c1==' ' && (c2=='P' || c2=='Q')) score+=-1.0;
    else if((c1=='P' || c1=='Q') && c2==' ') score+=-1.0;
    else if(c1=='P' && (c2=='A' || c2=='G' || c2=='R')) score+=1.0;
    else if(c1=='P' && (c2=='C' || c2=='U' || c2=='Y')) score+=0.5;
    else if((c1=='A' || c1=='G' || c1=='R') && c2=='P') score+=1.0;
    else if((c1=='C' || c1=='U' || c1=='Y') && c2=='P') score+=0.5;
    else if(c1=='Q' && (c2=='A' || c2=='G' || c2=='R')) score+=0.5;
    else if(c1=='Q' && (c2=='C' || c2=='U' || c2=='Y')) score+=1.0;
    else if((c1=='A' || c1=='G' || c1=='R') && c2=='Q') score+=0.5;
    else if((c1=='C' || c1=='U' || c1=='Y') && c2=='Q') score+=1.0;
    else {
      string str="match_RNA_motif.cpp: Unknown sequence symbol --> "+to_string(c1)+"/"+to_string(c2);
      error->one(str); 
    } 
  }

  score=score*1.0/seq1.size();
  if(score<0.0) score=0.0;
  
  return score;
}