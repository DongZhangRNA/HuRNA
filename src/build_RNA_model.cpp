#include "basic.h"
#include "initialization.h"
#include "error.h"
#include "read_PDB.h"
#include "measure.h"
#include "sort_vector.h"
#include "find_RNA_template.h"
#include "assemble_RNA_motif.h"
#include "build_RNA_model.h"

using namespace CBS_NS;

BuildRNAModel::BuildRNAModel(Cobmms *cbs) : Pointers(cbs), RNAHelix(cbs), RNALoop(cbs) {
  
}

BuildRNAModel::BuildRNAModel(const BuildRNAModel &other) : Pointers(other.cbs), RNAHelix(other.cbs), RNALoop(other.cbs) {
  models=other.models;
  adjusted_motifs=other.adjusted_motifs;
  Vtemplate_information=other.Vtemplate_information;

  RNAHelix=other.RNAHelix;
  RNALoop=other.RNALoop;
}

BuildRNAModel& BuildRNAModel::operator=(const BuildRNAModel &other) {
  if(this != &other) {
    Pointers::operator=(other.cbs);
    models=other.models;
    adjusted_motifs=other.adjusted_motifs;
    Vtemplate_information=other.Vtemplate_information;

    RNAHelix=other.RNAHelix;
    RNALoop=other.RNALoop;
  }

  return *this;
}

BuildRNAModel::~BuildRNAModel() {
  models.clear();
  models.shrink_to_fit();
  adjusted_motifs.clear();
  adjusted_motifs.shrink_to_fit();
  Vtemplate_information.clear();
  Vtemplate_information.shrink_to_fit();

  RNAHelix.clear();
  RNALoop.clear();
}

double BuildRNAModel::build_RNA_model(const vector<RNAMotif> &Vmotifs, ReadRestraint &restr_one) {
  models.clear();
  adjusted_motifs.clear();
  Vtemplate_information.clear();

  BuildRNAModel::read_RNA_helix_loop();

  vector<RNAMotif> Tmotifs;
  if(restr_one.RNA_templates.size()>0) {
    for(size_t ii=0;ii<Vmotifs.size();ii++) {
      int check=0;
      for(size_t jj=0;jj<restr_one.RNA_templates.size();jj++) {
        int relation=Measure::vector_relation(restr_one.RNA_templates[jj].Vresnames,Vmotifs[ii].Vresnames);
        if(relation==Measure::equal || relation==Measure::contained) {
          check=1;
          break;
        }
      }
      if(check==0) Tmotifs.push_back(Vmotifs[ii]);
    }
  }
  else {
    Tmotifs=Vmotifs;
  }

  int flag_junction=0;
  vector<RNAMotif> VLone_BP; //for Build_RNA_motif_structure
  for(size_t ii=0;ii<Tmotifs.size();ii++) {
    string str="build_RNA_model.cpp: Building motif "+to_string(ii+1)+"/"+to_string(Tmotifs.size())+" : ";
    str+=Tmotifs[ii].get_RNA_motif_info();
    //error->log(str);
    string type=Tmotifs[ii].type;
    int flag=-1;
    if(type=="Helix" || type=="Tertiary_helix" || type=="THelix") {
      BuildRNAModel::build_Helix(Tmotifs[ii]);
      flag=0;
    }
    else if(type=="Hairpin" || type=="Hairpin_kiss")
      flag=BuildRNAModel::build_Hairpin(Tmotifs[ii]);
    else if(type=="Tail5" || type=="Tail3")
      flag=BuildRNAModel::build_Tail(Tmotifs[ii]);
    else if(type.find("way_junction")!=string::npos) {
      flag=BuildRNAModel::build_Junction(Tmotifs[ii]);
      if(flag<0) flag_junction=1; //Failed to find template
    }
    else if(type=="Loop" || type=="Loop_kiss_loop") {
      BuildRNAModel::build_Loop(Tmotifs[ii]);
      flag=0;
    }
    else if(type=="Unpaired_chain") {
      BuildRNAModel::build_Long_loop(Tmotifs[ii]);
      flag=0;
    }
    else if(type=="Lone_BP") {
      VLone_BP.push_back(Tmotifs[ii]);
      continue;
    }
    else {
      string str="build_RNA_model.cpp: Unknown motif type --> "+type;
      error->one(str);
    }

    int flag_check=0;
    for(size_t jj=0;jj<restr_one.RNA_templates.size();jj++) {
      int relation=Measure::vector_relation(restr_one.RNA_templates[jj].Vresnames,Tmotifs[ii].Vresnames);
      if(relation==Measure::intersected) {
        flag_check=1;
        break;
      }
    }

    if(type=="2-way_junction" && Tmotifs[ii].lengths[0]<=1 && Tmotifs[ii].lengths[1]<=1) { //Do not treat as rigid body
      flag=0;
    }

    if(flag==1 && flag_check==0) {
      if(type.find("way_junction_kiss")==string::npos) {
        string resids=Tmotifs[ii].ends[0]+"-"+Tmotifs[ii].ends[1];
        for(size_t jj=1;jj<Tmotifs[ii].ends.size()/2;jj++)
          resids+=", "+Tmotifs[ii].ends[2*jj]+"-"+Tmotifs[ii].ends[2*jj+1];
        restr_one.add_rigid_body(resids);
      }
      else {
        size_t nends=Tmotifs[ii].ends.size();
        string resids=Tmotifs[ii].ends[nends-2]+", "+Tmotifs[ii].ends[nends-1];
        for(int jj=(int)nends/2-2;jj>=0;jj--)
          resids+=", "+Tmotifs[ii].ends[2*jj]+"-"+Tmotifs[ii].ends[2*jj+1];
        restr_one.add_rigid_body(resids);
      }
    }
  }

  if(flag_junction && VLone_BP.size()>0) {  //Lone_BP to close loops when failed to find template for N-way_junction
    for(size_t jj=0;jj<VLone_BP.size();jj++) {
      BuildRNAModel::build_Helix(VLone_BP[jj]);
    }
  }

  /*
  for(size_t ii=0;ii<Vtemplate_information.size();ii++) {
    cout<<Vtemplate_information[ii][0]<<endl;
    for(size_t jj=1;jj<Vtemplate_information[ii].size();jj++)
      cout<<"  "<<Vtemplate_information[ii][jj]<<endl;
  }
  */
  
  //Limit the number of models through limiting the number of templates
  map<string,size_t> Mtype_index;
  map<size_t,size_t> Midentical_motifs;
  vector<double> Vlength_eff;
  for(size_t ii=0;ii<adjusted_motifs.size();ii++) {
    string type=adjusted_motifs[ii][0].type;
    vector<string> &Vseq=adjusted_motifs[ii][0].sequences;
    vector<int> &Vlen=adjusted_motifs[ii][0].lengths;
    string type_one=type;
    for(size_t jj=0;jj<Vseq.size();jj++)
      type_one+="&"+Vseq[jj];
    double len_one=Vlen[0]*1.0;
    for(size_t jj=1;jj<Vlen.size();jj++)
      len_one+=Vlen[jj];
    if(type.find("2-way_junction")!=string::npos) {
      len_one=len_one*1.5;
      if(abs(Vlen[1]-Vlen[0])<=1) 
        len_one=len_one*0.5;
    }
    else if(type.find("way_junction")!=string::npos && type.find("2-way_junction")==string::npos)
      len_one=(len_one+Vseq.size())*Vseq.size()*1.0;
    else if(type=="Hairpin" || type=="Hairpin_kiss")
      len_one=len_one*0.5;
    else if(type=="Tail5" || type=="Tail3")
      len_one=len_one*0.3;
    else if(type=="Loop" || type=="Loop_kiss_loop" || type=="Unpaired_chain") 
      len_one=0.0;
    else 
      len_one=0.0;
    
    if(adjusted_motifs[ii].size()==1)
      len_one=0.0;

    if(Mtype_index.count(type_one)==0) {
      Mtype_index.emplace(type_one,ii);
    }
    else {
      Midentical_motifs.emplace(ii,Mtype_index[type_one]);
      len_one=0.0;
    }
    Vlength_eff.push_back(len_one);
  }
  Mtype_index.clear();

  vector<int> Vtemplate_used;
  int Nlarge=0, Nmedium=0, Nsmall=0;
  for(size_t ii=0;ii<Vlength_eff.size();ii++) {
    Vtemplate_used.push_back(1);
    if(Vlength_eff[ii]>=9) Nlarge++;
    else if(Vlength_eff[ii]>=1.5 && Vlength_eff[ii]<9) Nmedium++;
    else if(Vlength_eff[ii]>0 && Vlength_eff[ii]<1.5) Nsmall++;
  }

  vector<size_t> Vsort=SortVector::sort_vector_descending(Vlength_eff);
  if(Nlarge>=4) {
    Vtemplate_used[Vsort[0]]=2;
    Vtemplate_used[Vsort[1]]=2;
    Vtemplate_used[Vsort[2]]=2;
    Vtemplate_used[Vsort[3]]=2;
  }
  else if(Nlarge==3) {
    Vtemplate_used[Vsort[0]]=3;
    Vtemplate_used[Vsort[1]]=3;
    Vtemplate_used[Vsort[2]]=2;
  }
  else if(Nlarge==2) {
    Vtemplate_used[Vsort[0]]=4;
    Vtemplate_used[Vsort[1]]=4;
  }
  else if(Nlarge==1) {
    Vtemplate_used[Vsort[0]]=init->RNA_template_nmax;
    if(Nmedium>=2) {
      Vtemplate_used[Vsort[0]]=5;
      Vtemplate_used[Vsort[1]]=2;
      Vtemplate_used[Vsort[2]]=2;
    }
    else if(Nmedium==1)
      Vtemplate_used[Vsort[1]]=2;
  }
  else if(Nlarge==0) {
    if(Nmedium>=4) {
      Vtemplate_used[Vsort[0]]=2;
      Vtemplate_used[Vsort[1]]=2;
      Vtemplate_used[Vsort[2]]=2;
      Vtemplate_used[Vsort[3]]=2;
    }
    else if(Nmedium==3) {
      Vtemplate_used[Vsort[0]]=3;
      Vtemplate_used[Vsort[1]]=3;
      Vtemplate_used[Vsort[2]]=2;
    }
    else if(Nmedium==2) {
      Vtemplate_used[Vsort[0]]=4;
      Vtemplate_used[Vsort[1]]=4;
    }
    else if(Nmedium==1) {
      Vtemplate_used[Vsort[0]]=init->RNA_template_nmax;
      if(Nsmall>=2) {
        Vtemplate_used[Vsort[0]]=5;
        Vtemplate_used[Vsort[1]]=2;
        Vtemplate_used[Vsort[2]]=2;
      }
      else if(Nsmall==1)
        Vtemplate_used[Vsort[1]]=2;
    }
    else if(Nmedium==0) {
      if(Nsmall>=4) {
        Vtemplate_used[Vsort[0]]=2;
        Vtemplate_used[Vsort[1]]=2;
        Vtemplate_used[Vsort[2]]=2;
        Vtemplate_used[Vsort[3]]=2;
      }
      else if(Nsmall==3) {
        Vtemplate_used[Vsort[0]]=3;
        Vtemplate_used[Vsort[1]]=3;
        Vtemplate_used[Vsort[2]]=2;
      }
      else if(Nsmall==2) {
        Vtemplate_used[Vsort[0]]=4;
        Vtemplate_used[Vsort[1]]=4;
      }
      else if(Nsmall==1) {
        Vtemplate_used[Vsort[0]]=init->RNA_template_nmax;
      }
    }
  }

  for(size_t ii=0;ii<Vtemplate_used.size();ii++) {
    if(Vtemplate_used[ii] > init->RNA_template_nmax)
      Vtemplate_used[ii]=init->RNA_template_nmax;
  }

  for(auto & it : Midentical_motifs) {
    Vtemplate_used[it.first]=Vtemplate_used[it.second];
  }

  for(size_t ii=0;ii<adjusted_motifs.size();ii++) {
    int nm=static_cast<int>(adjusted_motifs[ii].size());
    if(nm>Vtemplate_used[ii]) {
      Tmotifs.clear();
      for(int jj=0;jj<Vtemplate_used[ii];jj++)
        Tmotifs.push_back(adjusted_motifs[ii][jj]);
      adjusted_motifs[ii]=Tmotifs;
    }
  }

  for(size_t ii=0;ii<Vtemplate_information.size();ii++) {
    int nm=static_cast<int> (Vtemplate_information[ii].size());
    if(nm<=2) continue;
    else if(nm>Vtemplate_used[ii]+1) {
      nm=Vtemplate_used[ii]+1;
    }
    for(int jj=1;jj<nm;jj++)
      Vtemplate_information[ii][jj]+=" (used)";
  }

  /*
  int n_HP=0;
  int n_2W=0;
  int n_NW=0;
  int n_TL=0;
  for(size_t ii=0;ii<adjusted_motifs.size();ii++) {
    if(adjusted_motifs[ii].size()==1) continue;
    string type=adjusted_motifs[ii][0].type;
    if(type=="Tail5" || type=="Tail3") n_TL++;
    else if(type=="Hairpin" || type=="Hairpin_kiss") n_HP++;
    else if(type.find("2-way_junction")!=string::npos) n_2W++;
    else if(type.find("way_junction")!=string::npos && type.find("2-way_junction")==string::npos) n_NW++;
    else {
      string str="build_RNA_model.cpp: Unknown motif type --> "+type;
      error->one(str);
    }
  }

  size_t n_HP_temp_max,n_2W_temp_max,n_NW_temp_max,n_TL_temp_max;
  if(n_NW<=1)
    n_NW_temp_max=init->RNA_template_nmax;
  else if(n_NW==2)
    n_NW_temp_max=3;
  else if(n_NW==3 || n_NW==4)
    n_NW_temp_max=2;
  else
    n_NW_temp_max=1;

  if(n_NW>0) {
    n_HP_temp_max=1;
    n_2W_temp_max=1;
    n_TL_temp_max=1;
  }
  else {
    if(n_2W<=1)
      n_2W_temp_max=init->RNA_template_nmax;
    else if(n_2W==2)
      n_2W_temp_max=3;
    else if(n_2W==3 || n_2W==4)
      n_2W_temp_max=2;
    else
      n_2W_temp_max=1;
    
    if(n_2W>0) {
      n_HP_temp_max=1;
      n_TL_temp_max=1;
    } 
    else {
      if(n_HP<=1)
        n_HP_temp_max=init->RNA_template_nmax;
      else if(n_HP==2)
        n_HP_temp_max=3;
      else if(n_HP==3 || n_HP==4)
        n_HP_temp_max=2;
      else
        n_HP_temp_max=1;

      if(n_HP>0)
        n_TL_temp_max=1;
      else {
        if(n_TL<=1)
          n_TL_temp_max=init->RNA_template_nmax;
        else if(n_TL==2)
          n_TL_temp_max=3;
        else if(n_TL==3 || n_TL==4)
          n_TL_temp_max=2;
        else
          n_TL_temp_max=1;
      }
    } 
  }
  */

  //cout<<"Adjusting..."<<endl;
  
  /*
  for(size_t ii=0;ii<adjusted_motifs.size();ii++) {
    size_t nm=adjusted_motifs[ii].size();
    if(nm==1) continue;
    string type=adjusted_motifs[ii][0].type;
    if(type=="Tail5" || type=="Tail3") {
      if(nm>n_TL_temp_max)
        nm=n_TL_temp_max;
    }
    else if(type=="Hairpin" || type=="Hairpin_kiss") {
      if(nm>n_HP_temp_max)
        nm=n_HP_temp_max;
    }
    else if(type.find("2-way_junction")!=string::npos) {
      if(nm>n_2W_temp_max)
        nm=n_2W_temp_max;
    }
    else if(type.find("way_junction")!=string::npos && type.find("2-way_junction")==string::npos) {
      if(nm>n_NW_temp_max)
        nm=n_NW_temp_max;
    }
    Tmotifs.clear();
    for(size_t jj=0;jj<nm;jj++)
      Tmotifs.push_back(adjusted_motifs[ii][jj]);
    adjusted_motifs[ii]=Tmotifs;
  }
  */
  
  if(restr_one.RNA_templates.size()>0) {
    size_t nRNAtemp=restr_one.RNA_templates.size();
    vector<size_t> Vsize;
    for(size_t ii=0;ii<nRNAtemp;ii++) {
      size_t num=restr_one.RNA_templates[0].number_of_residues();
      Vsize.push_back(num);
    }
    Vsort.clear();
    Vsort=SortVector::sort_vector_ascending(Vsize);
    for(size_t ii=0;ii<Vsort.size();ii++) {
      size_t index_i=Vsort[ii];
      vector<RNAMotif> Vtmp;
      Vtmp.push_back(restr_one.RNA_templates[index_i]);
      adjusted_motifs.insert(adjusted_motifs.begin(),Vtmp);
      Vtemplate_used.insert(Vtemplate_used.begin(),1);
      vector<string> temp_info;
      temp_info.push_back(restr_one.RNA_templates[index_i].get_RNA_motif_info());
      Vtemplate_information.insert(Vtemplate_information.begin(),temp_info);
    }

    map<size_t,size_t> Mtmp=Midentical_motifs;
    Midentical_motifs.clear();
    for(auto & it : Mtmp) {
      size_t index_1=it.first+nRNAtemp;
      size_t index_2=it.second+nRNAtemp;
      Midentical_motifs.emplace(index_1,index_2);
    }
    Mtmp.clear();
    /*
    vector<string> Vsize;
    for(size_t ii=0;ii<restr_one.RNA_templates.size();ii++) {
      string str=to_string(ii)+" "+to_string(restr_one.RNA_templates[0].number_of_residues());
      Vsize.push_back(str);
    }
    if(Vsize.size()>1)
      SortVector::sort_vector_ascending(Vsize,1);
    for(size_t ii=0;ii<Vsize.size();ii++){
      vector<string> stoken=split_string(Vsize[ii]);
      int index=stoi(stoken[0]);
      vector<RNAMotif> Vtmp;
      Vtmp.push_back(restr_one.RNA_templates[index]);
      adjusted_motifs.insert(adjusted_motifs.begin(),Vtmp);
      vector<string> temp_info;
      temp_info.push_back(restr_one.RNA_templates[index].get_RNA_motif_info());
      Vtemplate_information.insert(Vtemplate_information.begin(),temp_info);
    }
    */
  }

  //cout<<"Combining..."<<endl;
    
  vector<vector<size_t> > Vindexs;
  vector<size_t> Vacc;
  for(size_t ii=0;ii<adjusted_motifs.size();ii++) {
    size_t n0=1;
    for(size_t jj=ii+1;jj<adjusted_motifs.size();jj++) {
      if(Midentical_motifs.count(jj)==0)
        n0*=adjusted_motifs[jj].size();
      else 
        n0*=1;
    }
    Vacc.push_back(n0);
  }
  size_t num=0;
  while(num<Vacc[0]*adjusted_motifs[0].size()) {
    vector<size_t> indexs_one;
    for(size_t jj=0;jj<adjusted_motifs.size();jj++) {
      size_t nm=adjusted_motifs[jj].size();
      if(Midentical_motifs.count(jj))
        nm=1;
      size_t index=((size_t)(num/Vacc[jj]))%nm;
      indexs_one.push_back(index);
    }
    for(auto & it : Midentical_motifs)
      indexs_one[it.first]=indexs_one[it.second];
    num++;
    Vindexs.push_back(indexs_one);
  }

  /*
  for(size_t ii=0;ii<Vtemplate_information.size();ii++) {
    if(Vtemplate_information[ii].size()<=2) continue;
    vector<string> stoken=split_string(Vtemplate_information[ii][0]);
    string temp_type=stoken[0];
    if(temp_type.find("Hairpin")!=string::npos || temp_type.find("Tail")!=string::npos || temp_type.find("way_junction")!=string::npos) {
      vector<size_t> Vused;
      for(size_t jj=0;jj<Vindexs.size();jj++) {
        size_t index_i=Vindexs[jj][ii]+1;
        if(count(Vused.begin(),Vused.end(),index_i)==0)
          Vused.push_back(index_i);
      }
      for(size_t kk=1;kk<Vtemplate_information[ii].size();kk++) {
        if(count(Vused.begin(),Vused.end(),kk))
          Vtemplate_information[ii][kk]+=" (used)";
      }
    }
  }
  */

  //cout<<"Assembling..."<<endl;

  double moved_distance_max=-1.0;
  for(size_t ii=0;ii<Vindexs.size();ii++) {
    vector<vector<RNAMotif> > Vmodel_parts;
    vector<size_t> Vassembled;
    double moved_distance_max_i=-1.0;
    for(size_t jj=0;jj<adjusted_motifs.size();jj++) {
      if(adjusted_motifs[jj][0].type=="THelix") { //exclude THelix motif at the first round
        continue; 
      }
      if(count(Vassembled.begin(),Vassembled.end(),jj)) continue;
      vector<RNAMotif> Vpart;
      Vpart.push_back(adjusted_motifs[jj][Vindexs[ii][0]]);
      Vassembled.push_back(jj);
      vector<string> Vresnames_assembled=adjusted_motifs[jj][Vindexs[ii][0]].Vresnames;
      int flag_cycle=1;
      while (flag_cycle) {
        int flag_check=1;
        while(flag_check) {
          flag_check=0;
          for(size_t kk=0;kk<adjusted_motifs.size();kk++) {
            if(adjusted_motifs[kk][0].type=="THelix") continue;  //exclude THelix motif at the first round
            if(count(Vassembled.begin(),Vassembled.end(),kk)) continue;
            vector<string> &Vresnames_k=adjusted_motifs[kk][Vindexs[ii][kk]].Vresnames;
            int relation=Measure::vector_relation(Vresnames_assembled,Vresnames_k);
            if(relation==Measure::intersected || relation==Measure::contained || relation==Measure::equal) {
              flag_check=1;

              //for(size_t mm=0;mm<Vpart.size();mm++) {
                //cout<<" Ref: "<<Vpart[mm].get_RNA_motif_info()<<endl;
              //}
              //cout<<" Tar: "<<adjusted_motifs[kk][Vindexs[ii][kk]].get_RNA_motif_info()<<endl;
            
              AssembleRNAMotif ARmotif_one(cbs);
              ARmotif_one.assemble_RNA_motif(Vpart,adjusted_motifs[kk][Vindexs[ii][kk]],1);
              double moved_distance=ARmotif_one.moved_distance;
              if(moved_distance>moved_distance_max_i)
                moved_distance_max_i=moved_distance;
              Vpart.push_back(ARmotif_one.motif_adjust);
              Vassembled.push_back(kk);
              for(auto res_one : Vresnames_k) {
                if(count(Vresnames_assembled.begin(),Vresnames_assembled.end(),res_one)==0)
                  Vresnames_assembled.push_back(res_one);
              }
              break;
            }
          }
        }

        size_t nm_before=Vpart.size();

        //place THelix motif
        for(size_t kk=0;kk<adjusted_motifs.size();kk++) {
          if(adjusted_motifs[kk][0].type!="THelix") continue;
          if(count(Vassembled.begin(),Vassembled.end(),kk)) continue;
          vector<string> &Vresnames_k=adjusted_motifs[kk][Vindexs[ii][kk]].Vresnames;
          int relation=Measure::vector_relation(Vresnames_assembled,Vresnames_k);
          if(relation==Measure::intersected || relation==Measure::contained || relation==Measure::equal) {
            flag_check=1;

            //for(size_t mm=0;mm<Vpart.size();mm++) {
              //cout<<" Ref: "<<Vpart[mm].get_RNA_motif_info()<<endl;
            //}
            //cout<<" Tar: "<<adjusted_motifs[kk][Vindexs[ii][kk]].get_RNA_motif_info()<<endl;
          
            AssembleRNAMotif ARmotif_one(cbs);
            ARmotif_one.assemble_RNA_motif(Vpart,adjusted_motifs[kk][Vindexs[ii][kk]],1);
            double moved_distance=ARmotif_one.moved_distance;
            if(moved_distance>moved_distance_max_i)
              moved_distance_max_i=moved_distance;
            Vpart.push_back(ARmotif_one.motif_adjust);
            Vassembled.push_back(kk);
            for(auto res_one : Vresnames_k) {
              if(count(Vresnames_assembled.begin(),Vresnames_assembled.end(),res_one)==0)
                Vresnames_assembled.push_back(res_one);
            }
          }
        }
        
        if(Vpart.size()==nm_before) flag_cycle=0;
        if(Vassembled.size()==adjusted_motifs.size()) flag_cycle=0;
      } 
     
      string str="build_RNA_model.cpp: Building model "+to_string(ii+1)+"/"+to_string(Vindexs.size())+" --> ";
      str+="Assembled motifs ("+to_string(adjusted_motifs.size())+") :";
      for(auto nn : Vassembled)
        str+=" "+to_string(nn);
      //error->log(str);
      Vmodel_parts.push_back(Vpart);
    }

    string str="build_RNA_model.cpp: Building model "+to_string(ii+1)+"/"+to_string(Vindexs.size())+" --> ";
    str+="Move max "+to_string(moved_distance_max_i)+" Angstrom in placing motif.";
    //error->log(str);
    
    if(moved_distance_max<moved_distance_max_i)
      moved_distance_max=moved_distance_max_i;

    vector<RNAMotif> model_one=Vmodel_parts[0];
    if(Vmodel_parts.size()>1) {
      Vmodel_parts.erase(Vmodel_parts.begin());
      AssembleRNAMotif ARmotif_one(cbs);
      ARmotif_one.place_RNA_motif(model_one,Vmodel_parts);
      for(size_t kk=0;kk<Vmodel_parts.size();kk++) {
        model_one.insert(model_one.end(),Vmodel_parts[kk].begin(),Vmodel_parts[kk].end());
      }
    }
    models.push_back(model_one);
  }

  return moved_distance_max;
}

void BuildRNAModel::clear() {
  models.clear();
  models.shrink_to_fit();
  adjusted_motifs.clear();
  adjusted_motifs.shrink_to_fit();
  Vtemplate_information.clear();
  Vtemplate_information.shrink_to_fit();
}

void BuildRNAModel::read_RNA_helix_loop() {
  RNAHelix.clear();
  RNALoop.clear();

  string fstr=init->INPUT_DIR+"RNA_Aform_helix.pdb";
  ReadPDB readpdb_one(cbs);
  readpdb_one.read_PDB(fstr);
  RNAHelix.set_RNA_molecule(readpdb_one.RNA_atoms[0]);
  string str="RNAHelix: "+to_string(RNAHelix.chains.size())+" chains.";
  //error->message(str);
  max_helix_length=RNAHelix.chains[0].size();

  fstr=init->INPUT_DIR+"RNA_long_loop.pdb";
  readpdb_one.read_PDB(fstr);
  RNALoop.set_RNA_molecule(readpdb_one.RNA_atoms[0]);
  str="RNALoop: "+to_string(RNALoop.chains.size())+" chains.";
  //error->message(str);
  max_loop_length=RNALoop.chains[0].size();
}

void BuildRNAModel::build_Helix(RNAMotif &motif_one) {
  vector<string> temp_info;
  temp_info.push_back(motif_one.get_RNA_motif_info());
  size_t index_begin=0;
  size_t length=motif_one.chains[0].size();
  while(index_begin<length) {
    if((length-index_begin)<=max_helix_length) {
      if(motif_one.type=="Helix" || motif_one.type=="THelix") {
        Chain hx1(cbs),hx2(cbs);
        hx1.chain_label=Chain::rna;
        hx2.chain_label=Chain::rna;
        string secstr1="";
        string secstr2="";
        for(size_t ii=0;ii<length-index_begin;ii++) {
          hx1.insert_residue(motif_one.chains[0].residues[ii+index_begin],Chain::end);
          hx2.insert_residue(motif_one.chains[1].residues[length-ii-index_begin-1],Chain::begin);
          secstr1=secstr1+to_string((char)motif_one.secstrs[0][ii+index_begin]);   //"(";
          secstr2=to_string((char)motif_one.secstrs[1][length-ii-index_begin-1])+secstr2;   //")";
        }
        Chain hx1_xyz(cbs),hx2_xyz(cbs);
        hx1_xyz.chain_label=Chain::rna;
        hx2_xyz.chain_label=Chain::rna;
        for(size_t ii=0;ii<length-index_begin;ii++) {
          hx1_xyz.insert_residue(RNAHelix.chains[0].residues[ii],Chain::end);
          hx2_xyz.insert_residue(RNAHelix.chains[1].residues[max_helix_length-ii-1],Chain::begin);
        }
        hx1.set_chain_xyz(hx1_xyz);
        hx2.set_chain_xyz(hx2_xyz);
        RNAMotif helix_one(cbs);
        helix_one.type=motif_one.type;
        helix_one.insert_RNA_segment(hx1,0);
        helix_one.insert_RNA_segment(hx2,0);
        helix_one.secstrs.push_back(secstr1);
        helix_one.secstrs.push_back(secstr2);
        vector<RNAMotif> Vhelix;
        Vhelix.push_back(helix_one);
        adjusted_motifs.push_back(Vhelix);
        temp_info.push_back(helix_one.get_RNA_motif_info());
      }
      else if(motif_one.type=="Tertiary_helix") {
        Chain hx1(cbs);
        hx1.chain_label=Chain::rna;
        string secstr1="";
        for(size_t ii=0;ii<length-index_begin;ii++) {
          hx1.insert_residue(motif_one.chains[0].residues[ii+index_begin],Chain::end);
          secstr1=motif_one.secstrs[0].substr(index_begin);
        }
        Chain hx1_xyz(cbs);
        hx1_xyz.chain_label=Chain::rna;
        for(size_t ii=0;ii<length-index_begin;ii++) {
          hx1_xyz.insert_residue(RNAHelix.chains[0].residues[ii],Chain::end);
        }
        hx1.set_chain_xyz(hx1_xyz);
        RNAMotif Thelix_one(cbs);
        Thelix_one.type=motif_one.type;
        Thelix_one.insert_RNA_segment(hx1,0);
        Thelix_one.secstrs.push_back(secstr1);
        vector<RNAMotif> Vhelix;
        Vhelix.push_back(Thelix_one);
        adjusted_motifs.push_back(Vhelix);
        temp_info.push_back(Thelix_one.get_RNA_motif_info());
      }
      else if(motif_one.type=="Lone_BP") {
        Chain hx1(cbs),hx2(cbs);
        hx1.chain_label=Chain::rna;
        hx2.chain_label=Chain::rna;
        string secstr1=motif_one.secstrs[0];
        string secstr2=motif_one.secstrs[1];
        hx1.insert_residue(motif_one.chains[0].residues[0],Chain::begin);
        hx2.insert_residue(motif_one.chains[1].residues[0],Chain::begin);
        Chain hx1_xyz(cbs),hx2_xyz(cbs);
        hx1_xyz.chain_label=Chain::rna;
        hx2_xyz.chain_label=Chain::rna;
        hx1_xyz.insert_residue(RNAHelix.chains[0].residues[0],Chain::begin);
        hx2_xyz.insert_residue(RNAHelix.chains[1].residues.back(),Chain::begin);
        hx1.set_chain_xyz(hx1_xyz);
        hx2.set_chain_xyz(hx2_xyz);
        RNAMotif Lone_BP_one(cbs);
        Lone_BP_one.type=motif_one.type;
        Lone_BP_one.insert_RNA_segment(hx1,0);
        Lone_BP_one.insert_RNA_segment(hx2,0);
        Lone_BP_one.secstrs.push_back(secstr1);
        Lone_BP_one.secstrs.push_back(secstr2);
        vector<RNAMotif> Vhelix;
        Vhelix.push_back(Lone_BP_one);
        adjusted_motifs.push_back(Vhelix);
        temp_info.push_back(Lone_BP_one.get_RNA_motif_info());
      }
      else {
        string str="build_RNA_model.cpp: Not helix --> "+motif_one.type;
        error->one(str);
      }
      index_begin=length;
    }
    else {
      if(motif_one.type=="Helix" || motif_one.type=="THelix") {
        Chain hx1(cbs),hx2(cbs);
        hx1.chain_label=Chain::rna;
        hx2.chain_label=Chain::rna;
        string secstr1="";
        string secstr2="";
        for(size_t ii=0;ii<max_helix_length;ii++) {
          hx1.insert_residue(motif_one.chains[0].residues[ii+index_begin],Chain::end);
          hx2.insert_residue(motif_one.chains[1].residues[length-ii-index_begin-1],Chain::begin);
          secstr1=secstr1+to_string((char)motif_one.secstrs[0][ii+index_begin]);   //"(";
          secstr2=to_string((char)motif_one.secstrs[1][length-ii-index_begin-1])+secstr2;   //")";
        }
        Chain hx1_xyz(cbs),hx2_xyz(cbs);
        hx1_xyz.chain_label=Chain::rna;
        hx2_xyz.chain_label=Chain::rna;
        for(size_t ii=0;ii<max_helix_length;ii++) {
          hx1_xyz.insert_residue(RNAHelix.chains[0].residues[ii],Chain::end);
          hx2_xyz.insert_residue(RNAHelix.chains[1].residues[max_helix_length-ii-1],Chain::begin);
        }
        hx1.set_chain_xyz(hx1_xyz);
        hx2.set_chain_xyz(hx2_xyz);
        RNAMotif helix_one(cbs);
        helix_one.type=motif_one.type;
        helix_one.insert_RNA_segment(hx1,0);
        helix_one.insert_RNA_segment(hx2,0);
        helix_one.secstrs.push_back(secstr1);
        helix_one.secstrs.push_back(secstr2);
        vector<RNAMotif> Vhelix;
        Vhelix.push_back(helix_one);
        adjusted_motifs.push_back(Vhelix);
        temp_info.push_back(helix_one.get_RNA_motif_info());
      }
      else if(motif_one.type=="Tertiary_helix") {
        Chain hx1(cbs);
        hx1.chain_label=Chain::rna;
        string secstr1="";
        for(size_t ii=0;ii<max_helix_length;ii++) {
          hx1.insert_residue(motif_one.chains[0].residues[ii+index_begin],Chain::end);
          secstr1=motif_one.secstrs[0].substr(index_begin,max_helix_length);
        }
        Chain hx1_xyz(cbs);
        hx1_xyz.chain_label=Chain::rna;
        for(size_t ii=0;ii<max_helix_length;ii++) {
          hx1_xyz.insert_residue(RNAHelix.chains[0].residues[ii],Chain::end);
        }
        hx1.set_chain_xyz(hx1_xyz);
        RNAMotif Thelix_one(cbs);
        Thelix_one.type=motif_one.type;
        Thelix_one.insert_RNA_segment(hx1,0);
        Thelix_one.secstrs.push_back(secstr1);
        vector<RNAMotif> Vhelix;
        Vhelix.push_back(Thelix_one);
        adjusted_motifs.push_back(Vhelix);
        temp_info.push_back(Thelix_one.get_RNA_motif_info());
      }
      else {
        string str="build_RNA_model.cpp: Not helix --> "+motif_one.type;
        error->one(str);
      }
      index_begin+=max_helix_length-1;
    }
  }
  
  Vtemplate_information.push_back(temp_info);
  
}

void BuildRNAModel::build_Long_loop(RNAMotif &motif_one) {
  if(motif_one.chains.size()>1) {
    string str="build_RNA_model.cpp: More than one chain in build_Long_loop --> "+to_string(motif_one.chains.size());
    error->one(str);
  } 

  vector<string> temp_info;
  temp_info.push_back(motif_one.get_RNA_motif_info());
  size_t index_begin=0;
  size_t length=motif_one.chains[0].size();
  while(index_begin<length) {
    if((length-index_begin)<=max_loop_length) {
      Chain loop1(cbs);
      loop1.chain_label=Chain::rna;
      string secstr1="";
      for(size_t ii=0;ii<length-index_begin;ii++) {
        loop1.insert_residue(motif_one.chains[0].residues[ii+index_begin],Chain::end);
        secstr1=motif_one.secstrs[0].substr(index_begin);
      }
      Chain loop_xyz(cbs);
      loop_xyz.chain_label=Chain::rna;
      size_t index_offset=static_cast<size_t>(rand_one->uniform()*(max_loop_length-(length-index_begin)));
      for(size_t ii=0;ii<length-index_begin;ii++) {
        loop_xyz.insert_residue(RNALoop.chains[0].residues[ii+index_offset],Chain::end);
      }
      loop1.set_chain_xyz(loop_xyz);
      RNAMotif Lloop_one(cbs);
      Lloop_one.type="Long_loop";
      Lloop_one.insert_RNA_segment(loop1);
      Lloop_one.secstrs.push_back(secstr1);
      vector<RNAMotif> Vloop;
      Vloop.push_back(Lloop_one);
      adjusted_motifs.push_back(Vloop);
      temp_info.push_back(Lloop_one.get_RNA_motif_info());

      index_begin=length;
    }
    else {
      Chain loop1(cbs);
      loop1.chain_label=Chain::rna;
      string secstr1="";
      for(size_t ii=0;ii<max_loop_length;ii++) {
        loop1.insert_residue(motif_one.chains[0].residues[ii+index_begin],Chain::end);
        secstr1=motif_one.secstrs[0].substr(index_begin,max_loop_length);
      }
      Chain loop_xyz(cbs);
      loop_xyz.chain_label=Chain::rna;
      for(size_t ii=0;ii<max_loop_length;ii++) {
        loop_xyz.insert_residue(RNALoop.chains[0].residues[ii],Chain::end);
      }
      loop1.set_chain_xyz(loop_xyz);
      RNAMotif Lloop_one(cbs);
      Lloop_one.type="Long_loop";
      Lloop_one.insert_RNA_segment(loop1);
      Lloop_one.secstrs.push_back(secstr1);
      vector<RNAMotif> Vloop;
      Vloop.push_back(Lloop_one);
      adjusted_motifs.push_back(Vloop);
      temp_info.push_back(Lloop_one.get_RNA_motif_info());

      index_begin+=max_loop_length-1;
    }
  }

  Vtemplate_information.push_back(temp_info);
  
} 

int BuildRNAModel::build_Hairpin(RNAMotif &motif_one) {
  if(motif_one.type!="Hairpin" && motif_one.type!="Hairpin_kiss") {
    string str="build_RNA_model.cpp: Not Hairpin and Hairpin_kiss --> "+motif_one.type;
    error->one(str);
  } 

  vector<string> temp_info;
  temp_info.push_back(motif_one.get_RNA_motif_info());

  FindRNATemplate FRtemplate_one(cbs);
  FRtemplate_one.find_template(motif_one);
  if(FRtemplate_one.Vtemplates.size()>0) {
    if(FRtemplate_one.Vtemplate_scores[0].match_score==0) {
      vector<RNAMotif> Vmotif;
      for(size_t ii=0;ii<FRtemplate_one.Vtemplates.size();ii++) {
        RNAMotif motif_xyz=motif_one;
        motif_xyz.set_RNA_motif_xyz(FRtemplate_one.Vtemplates[ii]);
        Vmotif.push_back(motif_xyz);
        temp_info.push_back(FRtemplate_one.Vtemplate_scores[ii].get_match_info());
      }
      adjusted_motifs.push_back(Vmotif);
      Vtemplate_information.push_back(temp_info);

      if(FRtemplate_one.Vtemplate_scores[0].sequence_similarity>=init->RNA_sequence_similarity_high) { 
        size_t nlast=adjusted_motifs.size()-1;
        for(size_t jj=0;jj<adjusted_motifs[nlast].size();jj++) {
          adjusted_motifs[nlast][jj].flag_rigid=1;
        }
        return 1;
      }
      else
        return 0;
    }
    else {
      size_t length=motif_one.chains[0].size();
      int diff_length=(motif_one.lengths[0]-FRtemplate_one.Vtemplates[0].lengths[0])/2;
      Chain hx1(cbs),hx2(cbs),hp_loop(cbs);
      hx1.chain_label=Chain::rna;
      hx2.chain_label=Chain::rna;
      hp_loop.chain_label=Chain::rna;
      for(int ii=0;ii<=diff_length;ii++) {
        hx1.insert_residue(motif_one.chains[0].residues[ii],Chain::end);
        hx2.insert_residue(motif_one.chains[0].residues[length-ii-1],Chain::begin);
      }
      for(int ii=diff_length;ii<(int)length-diff_length;ii++)
        hp_loop.insert_residue(motif_one.chains[0].residues[ii],Chain::end);
      string secstr1=motif_one.secstrs[0].substr(0,diff_length+1);
      string secstr2=motif_one.secstrs[0].substr(diff_length,length-2*diff_length);
      string secstr3=motif_one.secstrs[0].substr(length-diff_length-1);
      RNAMotif hp_xyz(cbs);
      hp_xyz.type="Hairpin";
      hp_xyz.insert_RNA_segment(hp_loop,-2);
      hp_xyz.secstrs.push_back(secstr2);
      hp_xyz.set_RNA_motif_xyz(FRtemplate_one.Vtemplates[0]);
      vector<RNAMotif> Vmotif;
      Vmotif.push_back(hp_xyz);
      adjusted_motifs.push_back(Vmotif);
      temp_info.push_back(FRtemplate_one.Vtemplate_scores[0].get_match_info());
      RNAMotif helix_one(cbs);
      helix_one.type="Helix";
      helix_one.insert_RNA_segment(hx1);
      helix_one.insert_RNA_segment(hx2);
      helix_one.secstrs.push_back(secstr1);
      helix_one.secstrs.push_back(secstr3);
      temp_info.push_back(helix_one.get_RNA_motif_info());
      Vtemplate_information.push_back(temp_info);
      BuildRNAModel::build_Helix(helix_one);

      return 0;
    }
  }
  else {
    BuildRNAModel::build_Long_loop(motif_one);

    return 0;
  }
}

int BuildRNAModel::build_Tail(RNAMotif &motif_one) {
  if(motif_one.type!="Tail5" && motif_one.type!="Tail3") {
    string str="build_RNA_model.cpp: Not Tail --> "+motif_one.type;
    error->one(str);
  } 

  vector<string> temp_info;
  temp_info.push_back(motif_one.get_RNA_motif_info());

  FindRNATemplate FRtemplate_one(cbs);
  FRtemplate_one.find_template(motif_one);
  if(FRtemplate_one.Vtemplates.size()>0) {
    if(FRtemplate_one.Vtemplate_scores[0].match_score==0) {
      vector<RNAMotif> Vmotif;
      for(size_t ii=0;ii<FRtemplate_one.Vtemplates.size();ii++) {
        RNAMotif motif_xyz=motif_one;
        motif_xyz.set_RNA_motif_xyz(FRtemplate_one.Vtemplates[ii]);
        Vmotif.push_back(motif_xyz);
        temp_info.push_back(FRtemplate_one.Vtemplate_scores[ii].get_match_info());
      }
      adjusted_motifs.push_back(Vmotif);
      Vtemplate_information.push_back(temp_info);

      if(FRtemplate_one.Vtemplate_scores[0].sequence_similarity>=init->RNA_sequence_similarity_high) {
        size_t nlast=adjusted_motifs.size()-1;
        for(size_t jj=0;jj<adjusted_motifs[nlast].size();jj++) {
          adjusted_motifs[nlast][jj].flag_rigid=1;
        }
        return 1;
      }
      else
        return 0;
    }
    else {
      size_t length=FRtemplate_one.Vtemplates[0].chains[0].size();
      Chain loop_one(cbs);
      loop_one.chain_label=Chain::rna;
      for(size_t ii=0;ii<motif_one.chains[0].size();ii++) {
        if(motif_one.type=="Tail5")
          loop_one.insert_residue(FRtemplate_one.Vtemplates[0].chains[0].residues[length-ii-1],Chain::begin);
        else if(motif_one.type=="Tail3")
          loop_one.insert_residue(FRtemplate_one.Vtemplates[0].chains[0].residues[ii],Chain::end);
      }
      RNAMotif motif_xyz=motif_one;
      motif_xyz.chains[0].set_chain_xyz(loop_one);
      vector<RNAMotif> Vmotif;
      Vmotif.push_back(motif_xyz);
      adjusted_motifs.push_back(Vmotif);
      temp_info.push_back(FRtemplate_one.Vtemplate_scores[0].get_match_info());
      Vtemplate_information.push_back(temp_info);

      return 0;
    }
  }
  else {
    BuildRNAModel::build_Long_loop(motif_one);

    return 0;
  }
}

int BuildRNAModel::build_Junction(RNAMotif &motif_one) {
  if(motif_one.type.find("way_junction")==string::npos) {
    string str="build_RNA_model.cpp: Not N-way_junction --> "+motif_one.type;
    error->one(str);
  } 

  vector<string> temp_info;
  temp_info.push_back(motif_one.get_RNA_motif_info());

  size_t Nw=motif_one.lengths.size();
  FindRNATemplate FRtemplate_one(cbs);
  FRtemplate_one.find_template(motif_one);
  if(FRtemplate_one.Vtemplates.size()>0) {
    if(FRtemplate_one.Vtemplate_scores[0].match_score==0) {
      vector<RNAMotif> Vmotif;
      for(size_t ii=0;ii<FRtemplate_one.Vtemplates.size();ii++) {
        RNAMotif motif_xyz=motif_one;
        motif_xyz.set_RNA_motif_xyz(FRtemplate_one.Vtemplates[ii]);
        Vmotif.push_back(motif_xyz);
        temp_info.push_back(FRtemplate_one.Vtemplate_scores[ii].get_match_info());
      }
      adjusted_motifs.push_back(Vmotif);
      Vtemplate_information.push_back(temp_info);

      if(FRtemplate_one.Vtemplate_scores[0].sequence_similarity>=init->RNA_sequence_similarity_high) {
        size_t nlast=adjusted_motifs.size()-1;
        for(size_t jj=0;jj<adjusted_motifs[nlast].size();jj++) {
          adjusted_motifs[nlast][jj].flag_rigid=1;
        }
        return 1;
      }
      else
        return 0;
    }
    else {
      string match_code=FRtemplate_one.Vtemplate_scores[0].match_code;
      if(Nw!=match_code.size()) {
        string str="build_RNA_model.cpp: Number of chains are not matched --> "+to_string(Nw)+" : "+to_string(match_code.size());
        error->one(str);
      }
        
      RNAMotif motif_xyz(cbs);
      motif_xyz.type=motif_one.type+"_p";
      vector<RNAMotif> Vloops;
      for(size_t ii=0;ii<Nw;ii++) {
        if(match_code[ii]=='1') {
          Chain loop_i=motif_one.chains[ii];
          loop_i.set_chain_xyz(FRtemplate_one.Vtemplates[0].chains[ii]);
          motif_xyz.insert_RNA_segment(loop_i,-2);
          motif_xyz.secstrs.push_back(motif_one.secstrs[ii]);
        }
        else {
          RNAMotif loop_one(cbs);
          loop_one.type="Loop";
          loop_one.insert_RNA_segment(motif_one.chains[ii],-2);
          loop_one.secstrs.push_back(motif_one.secstrs[ii]);
          Vloops.push_back(loop_one);
        }
      }
      vector<RNAMotif> Vmotif;
      Vmotif.push_back(motif_xyz);
      adjusted_motifs.push_back(Vmotif);
      temp_info.push_back(FRtemplate_one.Vtemplate_scores[0].get_match_info());
      for(size_t ii=0;ii<Vloops.size();ii++)
        temp_info.push_back(Vloops[ii].get_RNA_motif_info());
      Vtemplate_information.push_back(temp_info);
      
      for(size_t ii=0;ii<Vloops.size();ii++)
        BuildRNAModel::build_Loop(Vloops[ii]);

      return -1;
    }
  }
  else {
    vector<RNAMotif> Vloops;
    for(size_t ii=0;ii<Nw;ii++) {
      RNAMotif loop_one(cbs);
      loop_one.type="Loop";
      loop_one.insert_RNA_segment(motif_one.chains[ii],-2);
      loop_one.secstrs.push_back(motif_one.secstrs[ii]);
      Vloops.push_back(loop_one);
      temp_info.push_back(loop_one.get_RNA_motif_info());
    }

    Vtemplate_information.push_back(temp_info);

    for(size_t ii=0;ii<Nw;ii++)
      BuildRNAModel::build_Loop(Vloops[ii]);

    return -1;
  }
}

void BuildRNAModel::build_Loop(RNAMotif &motif_one) {
  if(motif_one.type!="Loop" && motif_one.type!="Loop_kiss_loop") {
    string str="build_RNA_model.cpp: Not Loop and Loop_kiss_loop --> "+motif_one.type;
    error->one(str);
  } 

  vector<string> temp_info;
  temp_info.push_back(motif_one.get_RNA_motif_info());

  FindRNATemplate FRtemplate_one(cbs);
  FRtemplate_one.find_template(motif_one);
  if(FRtemplate_one.Vtemplates.size()>0) {
    vector<RNAMotif> Vmotif;
    for(size_t ii=0;ii<1;ii++) { //FRtemplate_one.Vtemplates.size()
      RNAMotif motif_xyz=motif_one;
      motif_xyz.set_RNA_motif_xyz(FRtemplate_one.Vtemplates[ii]);
      Vmotif.push_back(motif_xyz);
      temp_info.push_back(FRtemplate_one.Vtemplate_scores[ii].get_match_info());
    }
    adjusted_motifs.push_back(Vmotif);
    Vtemplate_information.push_back(temp_info);

    //if(FRtemplate_one.Vtemplate_scores[0].sequence_similarity>=init->RNA_sequence_similarity_high) 
      //return 1;
    //else
      //return 0;
  }
  else {
    BuildRNAModel::build_Long_loop(motif_one);

    //return 0;
  }
}

ostream& operator <<(ostream &output, BuildRNAModel &BRmodel_one) {
  output<<"List of available templates information for "<<BRmodel_one.Vtemplate_information.size()<<" motifs:"<<endl;
  for(size_t ii=0;ii<BRmodel_one.Vtemplate_information.size();ii++) {
    output<<BRmodel_one.Vtemplate_information[ii][0]<<endl;
    for(size_t jj=1;jj<BRmodel_one.Vtemplate_information[ii].size();jj++) 
      output<<"  "<<BRmodel_one.Vtemplate_information[ii][jj]<<endl;
    output<<endl;
  }

  return output;
}