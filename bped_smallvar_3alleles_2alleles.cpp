/*
 
copyright G. Nuel (2014-2022)

This program is distributed in the hope that it will be useful, but without any warranty. It is free to use for academic research, but not for commercial or medical purposes. For any question about this code, please contact the author (nuel@math.cnrs.fr).

build and test with:
 
g++ -std=c++0x -O3 bped_smallvar_3alleles_2alleles.cpp -o bped3alleles2alleles
./bped3alleles2alleles examples/test_3alleles2alleles.ped examples/test_3alleles2alleles_ev.txt

./bped3alleles2alleles examples/data-geno-Asc_fam9.ped

dot bnpedigree.dot -Tpdf > bnpedigree.pdf && open bnpedigree.pdf
 
dot jtpedigree.dot -Tpdf > jtpedigree.pdf && open jtpedigree.pdf
dot jtpedigree.dot -Grankdir=LR -Tpdf > jtpedigree.pdf && open jtpedigree.pdf
 
*/


#define SIZE 1

#define AFREQ1 0.084
#define AFREQ2 0.779
#define AFREQ3 0.137

#define BFREQ1 0.995
#define BFREQ2 0.005

#include "bnlib2.h"

variable alleleA(string name) {
    variable res(name,{"1","2","3"});
    return(res);
}
variable alleleB(string name) {
    variable res(name,{"1","2"});
    return(res);
}
variable genotypeA(string name) {
    variable res(name,{"11","21","31","12","22","32","13","23","33"});
    return(res);
}
variable genotypeB(string name) {
    variable res(name,{"11","21","12","22"});
    return(res);
}
variable genotype(string name) {
    variable res(name,{
        "1111","2111","3111","1211","2211","3211","1311","2311","3311",
        "1121","2121","3121","1221","2221","3221","1321","2321","3321",
        "1112","2112","3112","1212","2212","3212","1312","2312","3312",
        "1122","2122","3122","1222","2222","3222","1322","2322","3322"
        });
    return(res);
}

variable selector(string name) {
    variable res(name,{"pat","mat"});
    return(res);
}

// allele with maf
potential Kfounder_alleleA(string al,vector<double> freq={AFREQ1,AFREQ2,AFREQ3}) {

    potential res({alleleA(al)});
    map<string,string> conf;
    
    res[conf={{al,"1"}}]=freq[0];
    res[conf={{al,"2"}}]=freq[1];
    res[conf={{al,"3"}}]=freq[2];
    
    return res;
}
potential Kfounder_alleleB(string al,vector<double> freq={BFREQ1,BFREQ2}) {

    potential res({alleleB(al)});
    map<string,string> conf;
    
    res[conf={{al,"1"}}]=freq[0];
    res[conf={{al,"2"}}]=freq[1];
    
    return res;
}

// paternal/maternal selector for offspring meiosis at the marker
potential Koffspring_selector(string sel) {
    potential res({selector(sel)});
    map<string,string> conf;

    res[conf={{sel,"pat"}}]=0.5;
    res[conf={{sel,"mat"}}]=0.5;

    return(res);
}


// meiosis with selector for the disease allele
potential Koffspring_alleleA(string al,string pat,string mat,string sel) {
    
    variable v1=alleleA(al);
    variable v2=alleleA(pat);
    variable v3=alleleA(mat);
    variable v4=selector(sel);
    
    potential res({v1,v2,v3,v4});
    map<string,string> conf;
    
    for (auto i1 : v1.val()) for (auto i2 : v2.val()) for (auto i3 : v3.val()) for (auto i4 : v4.val()) {
            double val=0.;
            if ((i4)=="pat") {
                if (i1==i2)
                    val=1.;
            } else {
                if (i1==i3)
                    val=1.;
                    }
            conf[al]=i1; conf[pat]=i2; conf[mat]=i3; conf[sel]=i4;
            res[conf]=val;
        }
    
    return res;
}
potential Koffspring_alleleB(string al,string pat,string mat,string sel) {
    
    variable v1=alleleB(al);
    variable v2=alleleB(pat);
    variable v3=alleleB(mat);
    variable v4=selector(sel);
    
    potential res({v1,v2,v3,v4});
    map<string,string> conf;
    
    for (auto i1 : v1.val()) for (auto i2 : v2.val()) for (auto i3 : v3.val()) for (auto i4 : v4.val()) {
            double val=0.;
            if ((i4)=="pat") {
                if (i1==i2)
                    val=1.;
            } else {
                if (i1==i3)
                    val=1.;
                    }
            conf[al]=i1; conf[pat]=i2; conf[mat]=i3; conf[sel]=i4;
            res[conf]=val;
        }
    
    return res;
}

// deterministic genotype
potential KgenotypeA(string geno,string pat,string mat) {
    
    variable v1=genotypeA(geno);
    variable v2=alleleA(pat);
    variable v3=alleleA(mat);
    
    potential res({v1,v2,v3});
    map<string,string> conf;
    
    // initialize everything to 0.0
    for (auto i1 : v1.val()) for (auto i2 : v2.val()) for (auto i3 : v3.val()) {
        //cout<<"geno="<<i1<<" pat="<<i2<<" mat="<<i3<<endl;
        conf[geno]=i1; conf[pat]=i2; conf[mat]=i3;
        res[conf]=0.;
    }
    // and put the 1.0 where it should
   for (auto i2 : v2.val()) for (auto i3 : v3.val()) {
        //cout<<"geno="<<i1<<" pat="<<i2<<" mat="<<i3<<endl;
        conf[geno]=i2+i3; conf[pat]=i2; conf[mat]=i3;
        res[conf]=1.;
    }
    
    return res;
}
potential KgenotypeB(string geno,string pat,string mat) {
    
    variable v1=genotypeB(geno);
    variable v2=alleleB(pat);
    variable v3=alleleB(mat);
    
    potential res({v1,v2,v3});
    map<string,string> conf;
    
    // initialize everything to 0.0
    for (auto i1 : v1.val()) for (auto i2 : v2.val()) for (auto i3 : v3.val()) {
        //cout<<"geno="<<i1<<" pat="<<i2<<" mat="<<i3<<endl;
        conf[geno]=i1; conf[pat]=i2; conf[mat]=i3;
        res[conf]=0.;
    }
    // and put the 1.0 where it should
   for (auto i2 : v2.val()) for (auto i3 : v3.val()) {
        //cout<<"geno="<<i1<<" pat="<<i2<<" mat="<<i3<<endl;
        conf[geno]=i2+i3; conf[pat]=i2; conf[mat]=i3;
        res[conf]=1.;
    }
    
    return res;
}
potential Kgenotype(string geno,string genoA,string genoB) {

    variable v1=genotype(geno);
    variable v2=genotypeA(genoA);
    variable v3=genotypeB(genoB);
    
    potential res({v1,v2,v3});
    map<string,string> conf;

    // initialize everything to 0.0
    for (auto i1 : v1.val()) for (auto i2 : v2.val()) for (auto i3 : v3.val()) {
        //cout<<"geno="<<i1<<" pat="<<i2<<" mat="<<i3<<endl;
        conf[geno]=i1; conf[genoA]=i2; conf[genoB]=i3;
        res[conf]=0.;
    }
    // and put the 1.0 where it should
    for (auto i2 : v2.val()) for (auto i3 : v3.val()) {
        //cout<<"geno="<<i1<<" pat="<<i2<<" mat="<<i3<<endl;
        conf[geno]=i2+i3; conf[genoA]=i2; conf[genoB]=i3;
        res[conf]=1.;
    }

    return res;
}

int main (int argc, char *argv[]) {
    
    vector<double> Afreq={AFREQ1,AFREQ2,AFREQ3};
    vector<double> Bfreq={BFREQ1,BFREQ2};

    // load variables and potentials
    map<string,variable> V;
    vector<string> id;
    map<string,potential> K;
    
    // read freq
    if (argc>=8) {
        string sAfreq1(argv[3]);
        string sAfreq2(argv[4]);
        string sAfreq3(argv[5]);
        string sBfreq1(argv[6]);
        string sBfreq2(argv[7]);
        if (sAfreq1!="") Afreq[0]=stod(sAfreq1);
        if (sAfreq2!="") Afreq[1]=stod(sAfreq2);
        if (sAfreq3!="") Afreq[2]=stod(sAfreq3);
        if (sBfreq1!="") Bfreq[0]=stod(sBfreq1);
        if (sBfreq2!="") Bfreq[1]=stod(sBfreq2);
    }
    
    // read .ped input
    {
        string famid,indid,patid,matid,gender,twin,pgm1,radin,rh,empty;
        ifstream myfile (argv[1]);
        if (myfile.is_open()) {
            while ( !myfile.eof() ) {
                getline (myfile,famid,'\t');
                getline (myfile,indid,'\t');
                getline (myfile,patid,'\t');
                getline (myfile,matid,'\n');
                if (famid!="") {
                    indid=famid+":"+indid;
                    id.push_back(indid);
                    //cout<<"processing "<<indid<<endl;
                    // add alleles
                    V["AXpat_"+indid]=alleleA("AXpat_"+indid);
                    V["AXmat_"+indid]=alleleA("AXmat_"+indid);
                    V["BXpat_"+indid]=alleleB("BXpat_"+indid);
                    V["BXmat_"+indid]=alleleB("BXmat_"+indid);
                    if ((patid!="0")&(matid!="0")) {
                        // offspring
                        patid=famid+":"+patid;
                        matid=famid+":"+matid;
                        //cout<<"pat="<<patid<<" mat="<<matid<<endl;
                        // add selectors
                        V["SAXpat_"+indid]=selector("SAXpat_"+indid);
                        V["SAXmat_"+indid]=selector("SAXmat_"+indid);
                        V["SBXpat_"+indid]=selector("SBXpat_"+indid);
                        V["SBXmat_"+indid]=selector("SBXmat_"+indid);
                        // and their potentials
                        K["SAXpat_"+indid]=Koffspring_selector("SAXpat_"+indid);
                        K["SAXmat_"+indid]=Koffspring_selector("SAXmat_"+indid);
                        K["AXpat_"+indid]=Koffspring_alleleA("AXpat_"+indid,"AXpat_"+patid,"AXmat_"+patid,"SAXpat_"+indid);
                        K["AXmat_"+indid]=Koffspring_alleleA("AXmat_"+indid,"AXpat_"+matid,"AXmat_"+matid,"SAXmat_"+indid);
                        K["SBXpat_"+indid]=Koffspring_selector("SBXpat_"+indid);
                        K["SBXmat_"+indid]=Koffspring_selector("SBXmat_"+indid);
                        K["BXpat_"+indid]=Koffspring_alleleB("BXpat_"+indid,"BXpat_"+patid,"BXmat_"+patid,"SBXpat_"+indid);
                        K["BXmat_"+indid]=Koffspring_alleleB("BXmat_"+indid,"BXpat_"+matid,"BXmat_"+matid,"SBXmat_"+indid);
                    } else {
                        // founder
                        K["AXpat_"+indid]=Kfounder_alleleA("AXpat_"+indid,Afreq);
                        K["AXmat_"+indid]=Kfounder_alleleA("AXmat_"+indid,Afreq);
                        K["BXpat_"+indid]=Kfounder_alleleB("BXpat_"+indid,Bfreq);
                        K["BXmat_"+indid]=Kfounder_alleleB("BXmat_"+indid,Bfreq);
                    }
                    // add genotype
                    V["AX_"+indid]=genotypeA("AX_"+indid);
                    K["AX_"+indid]=KgenotypeA("AX_"+indid,"AXpat_"+indid,"AXmat_"+indid);
                    V["BX_"+indid]=genotypeB("BX_"+indid);
                    K["BX_"+indid]=KgenotypeB("BX_"+indid,"BXpat_"+indid,"BXmat_"+indid);
                    V["X_"+indid]=genotype("X_"+indid);
                    K["X_"+indid]=Kgenotype("X_"+indid,"AX_"+indid,"BX_"+indid);
                }
            }
            myfile.close();
        }
        else cout << "Unable to open the pedigree file"<<endl;
    }
    
    // read evidence file
    if (argc>=3) {
        string famid,indid;
        // B=11
        string ev1111,ev2111,ev3111,ev1211,ev2211,ev3211,ev1311,ev2311,ev3311;
        // B=21
        string ev1121,ev2121,ev3121,ev1221,ev2221,ev3221,ev1321,ev2321,ev3321;
        // B=12
        string ev1112,ev2112,ev3112,ev1212,ev2212,ev3212,ev1312,ev2312,ev3312;
        // B=22
        string ev1122,ev2122,ev3122,ev1222,ev2222,ev3222,ev1322,ev2322,ev3322;
        ifstream myfile (argv[2]);
        if (myfile.is_open()) {
            while ( !myfile.eof() ) {
                getline (myfile,famid,'\t');
                getline (myfile,indid,'\t');
                // B=11
                getline (myfile,ev1111,'\t');
                getline (myfile,ev2111,'\t');
                getline (myfile,ev3111,'\t');
                getline (myfile,ev1211,'\t');
                getline (myfile,ev2211,'\t');
                getline (myfile,ev3211,'\t');
                getline (myfile,ev1311,'\t');
                getline (myfile,ev2311,'\t');
                getline (myfile,ev3311,'\t');
                // B=21
                getline (myfile,ev1121,'\t');
                getline (myfile,ev2121,'\t');
                getline (myfile,ev3121,'\t');
                getline (myfile,ev1221,'\t');
                getline (myfile,ev2221,'\t');
                getline (myfile,ev3221,'\t');
                getline (myfile,ev1321,'\t');
                getline (myfile,ev2321,'\t');
                getline (myfile,ev3321,'\t');
                // B=12
                getline (myfile,ev1112,'\t');
                getline (myfile,ev2112,'\t');
                getline (myfile,ev3112,'\t');
                getline (myfile,ev1212,'\t');
                getline (myfile,ev2212,'\t');
                getline (myfile,ev3212,'\t');
                getline (myfile,ev1312,'\t');
                getline (myfile,ev2312,'\t');
                getline (myfile,ev3312,'\t');
                // B=22
                getline (myfile,ev1122,'\t');
                getline (myfile,ev2122,'\t');
                getline (myfile,ev3122,'\t');
                getline (myfile,ev1222,'\t');
                getline (myfile,ev2222,'\t');
                getline (myfile,ev3222,'\t');
                getline (myfile,ev1322,'\t');
                getline (myfile,ev2322,'\t');
                getline (myfile,ev3322,'\n');
                //cout<<"famid="<<famid<<" indid="<<indid<<" ev3311="<<ev3311<<" ev3322="<<ev3322<<endl;
                if (famid!="") {
                    indid=famid+":"+indid;
                    //cout<<"G_"<<indid<<endl;
                    //V["X_"+indid]=genotype("X_"+indid);
                    // evidence
                    potential ev({V["X_"+indid]});
                    map<string,string> conf;
                    // B=11
                    ev[conf={{"X_"+indid,"1111"}}]=stod(ev1111);
                    ev[conf={{"X_"+indid,"2111"}}]=stod(ev2111);
                    ev[conf={{"X_"+indid,"3111"}}]=stod(ev3111);
                    ev[conf={{"X_"+indid,"1211"}}]=stod(ev1211);
                    ev[conf={{"X_"+indid,"2211"}}]=stod(ev2211);
                    ev[conf={{"X_"+indid,"3211"}}]=stod(ev3211);
                    ev[conf={{"X_"+indid,"1311"}}]=stod(ev1311);
                    ev[conf={{"X_"+indid,"2311"}}]=stod(ev2311);
                    ev[conf={{"X_"+indid,"3311"}}]=stod(ev3311);
                    // B=21
                    ev[conf={{"X_"+indid,"1121"}}]=stod(ev1121);
                    ev[conf={{"X_"+indid,"2121"}}]=stod(ev2121);
                    ev[conf={{"X_"+indid,"3121"}}]=stod(ev3121);
                    ev[conf={{"X_"+indid,"1221"}}]=stod(ev1221);
                    ev[conf={{"X_"+indid,"2221"}}]=stod(ev2221);
                    ev[conf={{"X_"+indid,"3221"}}]=stod(ev3221);
                    ev[conf={{"X_"+indid,"1321"}}]=stod(ev1321);
                    ev[conf={{"X_"+indid,"2321"}}]=stod(ev2321);
                    ev[conf={{"X_"+indid,"3321"}}]=stod(ev3321);
                    // B=12
                    ev[conf={{"X_"+indid,"1112"}}]=stod(ev1112);
                    ev[conf={{"X_"+indid,"2112"}}]=stod(ev2112);
                    ev[conf={{"X_"+indid,"3112"}}]=stod(ev3112);
                    ev[conf={{"X_"+indid,"1212"}}]=stod(ev1212);
                    ev[conf={{"X_"+indid,"2212"}}]=stod(ev2212);
                    ev[conf={{"X_"+indid,"3212"}}]=stod(ev3212);
                    ev[conf={{"X_"+indid,"1312"}}]=stod(ev1312);
                    ev[conf={{"X_"+indid,"2312"}}]=stod(ev2312);
                    ev[conf={{"X_"+indid,"3312"}}]=stod(ev3312);
                    // B=22
                    ev[conf={{"X_"+indid,"1122"}}]=stod(ev1122);
                    ev[conf={{"X_"+indid,"2122"}}]=stod(ev2122);
                    ev[conf={{"X_"+indid,"3122"}}]=stod(ev3122);
                    ev[conf={{"X_"+indid,"1222"}}]=stod(ev1222);
                    ev[conf={{"X_"+indid,"2222"}}]=stod(ev2222);
                    ev[conf={{"X_"+indid,"3222"}}]=stod(ev3222);
                    ev[conf={{"X_"+indid,"1322"}}]=stod(ev1322);
                    ev[conf={{"X_"+indid,"2322"}}]=stod(ev2322);
                    ev[conf={{"X_"+indid,"3322"}}]=stod(ev3322);
                    // and the potential
                    //K["X_"+indid]=Kgenotype("X_"+indid,"Xpat_"+indid,"Xmat_"+indid)*ev;
                    K["X_"+indid]=K["X_"+indid]*ev;
                }
            }
            myfile.close();
        }
        else cout << "Unable to open the evidence file"<<endl;
    }

    
    // build BN
    bn mybn;
    for (map<string,variable>::iterator it=V.begin(); it!=V.end(); it++) {
        //cout<<"V["<<it->first<<"]"<<endl;
        mybn.add(it->second,K[it->first]);
    }
    // validate and display BN
    mybn.initialize(false);
    //mybn.bndot("bnpedigree.dot");
    
    
    // build and display JT
    mybn.jt(false);
    //mybn.jtdot("jtpedigree.dot");
    
    // perform belief propagation
    mybn.bp(false);
    
    // posterior distribution
    for (auto myid : id) {
        cout<<"X_"<<myid<<"\t";
        mybn.posterior(V["X_"+myid]).print_data(cout);
        cout<<endl;
    }
    
    //cout<<"loglik(freq="<<freq[0]<<","<<freq[1]<<","<<freq[2]<<")="<<mybn.pevidence()<<endl;
    
    return 0;
};


