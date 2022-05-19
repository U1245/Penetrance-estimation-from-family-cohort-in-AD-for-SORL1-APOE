/*
 
 bnlib - a C++ Bayesian Network Library
 
 copyright: G. Nuel (2014-2022).

This program is distributed in the hope that it will be useful, but without any warranty. It is free to use for academic research, but not for commercial or medical purposes. For any question about this code, please contact the author (nuel@math.cnrs.fr).
 
 ChangeLog:
 - 0.1: basic version with full propagation
 - 0.2: data in potentials now in log scale
 - 0.3: potential now contains valarray
 - 0.4: major code overhaul (individual variable class,
 potential now created with set of variables)
 - 0.5: major logscale revamp.
 
 Known bugs:
 - marginal of impossible potential is incorrect
 - test for relerr wrong when log(_pevidence)=0.0

 Manual:
 
 1) put the following in your file mybn.cpp (and adjust SIZE to fit your needs):
 
 #define SIZE 3
 #include "bnlib.h"
 
 2) build and run with:
 
 g++ -std=c++0x -O3 mybn.cpp -o mybn && ./mybn
 
 3) BN and JT with:
 
 dot bn.dot -Tpdf > bn.pdf && open bn.pdf
 
 dot jt.dot -Tpdf > jt.pdf && open jt.pdf
 
 */

#ifndef BNLIB_H
#define BNLIB_H


#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <set>
#include <climits>
#include <cmath>
#include <valarray>
using namespace std;

// basic output functions
ostream &operator<<(ostream &output,const valarray<double> &val) {
    for (int i=0; i<val.size(); i++)
        output<<val[i]<<" ";
    return output;
}
ostream &operator<<(ostream &output,const vector<int>& vect) {
    for (int i=0; i<vect.size(); i++)
        output<<vect[i]<<" ";
    return output;
}
ostream &operator<<(ostream &output,const set<int>& S) {
    for (set<int>::const_iterator it=S.begin(); it!=S.end(); it++)
        output<<*it<<" ";
    return output;
}


// very simple variable class with name and string values
// only finite discrete variables are supported
class variable {
    friend class potential;
    friend class bn;
    
    string _name;
    vector<string> _val;
    string _id;
    
    void update(string val) {
        _id+=val;
    }
    void update() {
        _id=_name;
        for (vector<string>::iterator it=_val.begin(); it!=_val.end(); it++)
            update(*it);
    }
    
    
public:
    
    variable() {
        _name="undef";
        update();
    }
    
    variable(string name,vector<string> val) {
        _name=name;
        _val=val;
        update();
    }
    variable(string name) {
        _name=name;
        update();
    }
    variable(const variable&ref,string newname) {
        _name=newname;
        _val=ref._val;
        update();
    }
    variable& rename(string name) {
        _name=name;
        update();
        return *this;
    }
    variable& add_val(string val) {
        _val.push_back(val);
        update(val);
        return *this;
    }
    
    inline void operator=(const variable &var)
    {
        _name=var._name;
        _val=var._val;
        _id=var._id;
    }
    
    inline string& id() {
        return _id;
    }
    inline string& name() {
        return _name;
    }
    inline string name() const {
        return _name;
    }
    inline int dim() const {
        return _val.size();
    }
    inline vector<string> val() const {
        return _val;
    }
    inline bool operator==(const variable& var) const {
        return(this->_id==var._id);
    }
    inline bool operator!=(const variable& var) const {
        return(this->_id!=var._id);
    }
    inline bool operator<(const variable& var) const {
        return(this->_id<var._id);
    }
    inline bool operator>(const variable& var) const {
        return(this->_id>var._id);
    }
    inline bool operator()(const variable& var1,const variable& var2 ) const {
        return(var1._id>var2._id);
    }
    friend ostream &operator<<( ostream &output,const variable &var) {
        output << var._name + " = ";
        if (var.dim()==0)
            output<<"NULL";
        for (vector<string>::const_iterator it=var._val.begin(); it!=var._val.end(); it++) {
            output<<*it;
            if ((it+1)!=var._val.end()) output<<" or ";
        }
        return output;
    }
    
};

// sophisticated potential class
// support log-scale
// logscale is now a valarray<double> with predefind SIZE
// uses valarray<double> with predefined SIZE
class potential {
    friend class bn;
    
    // general
    int _size;
    // variables
    vector<variable> _var;
    // index
    vector<int> _dim;
    vector<int> _dim_index;
    vector<int> _current_conf;
    int _current_index;
    // data
    vector<valarray<double> > _data;
    valarray<double> _sum;
    valarray<double> _logscale;
    
    // reset the index
    void reset() {
        _current_index=0;
        for (vector<int>::iterator it=_current_conf.begin(); it!=_current_conf.end(); it++)
            *it=0;
    }
    
    int increment_index(int rank) {
        if (rank==-1) return -1;
        if (rank>=_size) return -1;
        if (_size==0) return -1;
        if (rank>-1) {
            _current_conf[rank]++;
            if (_current_conf[rank]==_dim[rank]) {
                _current_conf[rank]=0;
                _current_index-=(_dim[rank]-1)*_dim_index[rank];
            } else {
                _current_index+=_dim_index[rank];
            }
        }
        return _current_conf[rank];
    }
    
    void print_current_conf(ostream &output) const {
        if (_size==0) {
            output<<"NULL -> "<<_sum<<endl;
            return;
        }
        //cout<<_current_index<<": ";
        for (vector<int>::const_iterator it=_current_conf.begin(); it!=_current_conf.end(); it++) {
            int i=distance(_current_conf.begin(),it);
            output<<_var[i].name()<<"="<<_var[i]._val[*it];
            output<<" ";
        }
        output<<"-> "<<_data[_current_index]*exp(_logscale);
        output<<endl;
    }
    
    // allocate workspace and compute _dim_index
    void initialize() {
        _size=_var.size();
        _dim_index.clear();
        int aux=1;
        for (vector<int>::iterator it=_dim.begin(); it!=_dim.end(); it++) {
            _dim_index.push_back(aux);
            aux*=*it;
        }
        valarray<double> tmp(0.,SIZE); 
        _data.resize(aux,tmp);
        if (_logscale.size()!=SIZE) _logscale.resize(SIZE,0.);
        _current_conf.resize(_size,0.0);
    }
    
    // return the potential value at the current index
    inline valarray<double> &get() {
        if (_size>0)
            return _data[_current_index];
        else
            return _sum;
    }
  
    
public:
    
    inline valarray<double> logscale() const {
        return _logscale;
    }
    
    // crude data output in the case where SIZE=1
    void print_data(ostream &output) {
        for (vector<valarray<double> >::iterator it=_data.begin(); it!=_data.end(); it++)
            output<<(*it)[0]*exp(_logscale)<<"\t";
    }
    
    // compute, store and return the potential sum
    valarray<double> & sum() {
        if (_size==0) return _sum;
        _sum.resize(SIZE,0.);
        for (vector<valarray<double> >::iterator it=_data.begin(); it!=_data.end(); it++)
                _sum+=*it;
        return _sum;
    }
    
    // rescale the potential such as _sum=1.
    void rescale() {
        // compute the sum
        sum();
        // update the logscale
        _logscale+=log(_sum);
        // rescale the potential
        for (vector<valarray<double> >::iterator it=_data.begin(); it!=_data.end(); it++)
            *it/=_sum;
        // update the sum
        _sum=1.;
    }
    
    // empty constructor
    potential() {
        _size=0;
        _sum.resize(SIZE,0.);
        _logscale.resize(SIZE,0.);
    }
    
    // build a potential from a set of variables
    potential(const set<variable>& var) {
        _size=0;
        if (var.size()>0) {
            int aux=1;
            for (set<variable>::const_iterator it=var.begin(); it!=var.end(); it++) {
                _var.push_back(*it);
                _dim.push_back(it->dim());
                _dim_index.push_back(aux);
                aux*=_dim[_size++];
            }
            initialize();
        }
        _current_index=0;
    }
    
    // potential output
    friend ostream &operator<<(ostream &output,potential pot) {
        if (pot._size==0)
            pot.print_current_conf(output);
        pot.reset();
        int rank=0;
        while (rank<pot._size) {
            if (pot.get().max()>0.)
                    pot.print_current_conf(output);
            rank=0;
            while (pot.increment_index(rank)==0) rank++;
        }
        return output;
    }
    
    // subscripting operator for potentials
    valarray<double>& operator[](map<string,string> conf) {
        if (_size==0) {
            cerr<<"operator[] error: empty potential"<<endl;
            exit(0);
        }
        // locate conf
        _current_index=0;
        for (vector<variable>::iterator it1=_var.begin(); it1!=_var.end(); it1++) {
            map<string,string>::iterator it2=conf.find(it1->name());
            if (it2==conf.end()) {
                cerr<<"operator[] error: '"<<it1->name()<<"' is not provided in the configuration"<<endl;
                exit(0);
            }
            vector<string>::iterator it3;
            for (it3=it1->_val.begin(); it3!=it1->_val.end(); it3++) {
                if (*it3==it2->second) break;
            }
            if (it3==it1->_val.end()) {
                cerr<<"operator[] error: '"<<it2->second<<"' is not a valid value for variable '"<<it2->first<<"'"<<endl;
                exit(0);
            }
            int i=distance(_var.begin(),it1);
            int j=distance(it1->_val.begin(),it3);
            _current_conf[i]=j;
            _current_index+=_dim_index[i]*j;
        }
        return(_data[_current_index]);
    }
    
    // build a new potential by marginalizing the current potential
    potential marginal(const set<variable>& downto) {
        // initialize result potential
        potential res;
        res._logscale=_logscale;
        if (downto.size()==0) {
            res._sum=sum();
            return res;
        }
        // copy downto set
        for (set<variable>::iterator it=downto.begin(); it!=downto.end(); it++)
            res._var.push_back(*it);
        res._size=res._var.size();
        if (res._size>_size) {
            cerr<<"call of marginal where 'downto' is not a valid subset"<<endl;
            return(potential());
        }
        // build rank_mapping such as rank_mapping[i]=-1 if variable i in potential is not in downto
        // and rank_mapping[i]=j if variable i is variable j in the new potential
        vector<int> rank_mapping(_size,-1);
        int pos=0;
        for (int i=0; i<_size; i++) {
            if (pos<res._var.size() && _var[i]==res._var[pos]) {
                rank_mapping[i]=pos++;
                res._dim.push_back(_dim[i]);
            }
        }
        // check consistency
        if (pos<res._size) {
            cerr<<"call of marginal where 'downto' is not a valid subset"<<endl;
            return(potential());
        }
        res._size=pos;
        
        // initialize potential workspace
        res.initialize();
        //res._logscale=_logscale;
        
        // main loop
        reset(); res.reset();
        int rank=0;
        while (rank<_size) {
            res._data[res._current_index]+=_data[_current_index];
            rank=0;
            res.increment_index(rank_mapping[rank]);
            while (increment_index(rank)==0) {
                rank++;
                if (rank<_size) res.increment_index(rank_mapping[rank]);
            }
        }
        
        // return result
        return(res);
    }
    
    // product of two potential
    // fixme for performance
    //potential operator*(potential & two) {
    potential operator*(potential two) {
        potential one=*this;
        // initialize result potential
        potential res;
        res._logscale=one._logscale+two._logscale;
        // particular case.
        if (one._size==0) {
            //cout<<"one"<<endl;
            res=two;
            res._logscale=one._logscale+two._logscale;
            //cout<<"res._data.size()="<<res._data.size()<<endl;
            for (vector<valarray<double> >::iterator it=res._data.begin(); it!=res._data.end(); it++)
                *it*=one._sum;
            res._sum*=one._sum;
            //cout<<"ok"<<endl;
            return res;
        }
        if (two._size==0) {
            res=one;
            res._logscale=one._logscale+two._logscale;
            for (vector<valarray<double> >::iterator it=res._data.begin(); it!=res._data.end(); it++)
                *it*=two._sum;
            res._sum*=two._sum;
            return res;
        }
        // particular case
        if ((one._size==0) && (two._size==0)) {
            //cout<<"special case :"<<endl;
            res._sum=one.get()*two.get();
            //cout<<res.get()<<" "<<one.get()<<" "<<two.get()<<endl;
            return res;
        }
        set<variable> aux;
        aux.insert(one._var.begin(),one._var.end());
        aux.insert(two._var.begin(),two._var.end());
        for (set<variable>::iterator it=aux.begin(); it!=aux.end(); it++) {
            res._var.push_back(*it);
            res._dim.push_back(it->dim());
        }
        res.initialize();
        //cout<<"res._size="<<res._size<<endl;
        
        // compute rank mappings
        vector<int> rank_mapping1(res._size,-1);
        vector<int> rank_mapping2(res._size,-1);
        {
            int i1=0,i2=0;
            for (int j=0; j<res._size; j++) {
                if (i1<one._var.size() && res._var[j]==one._var[i1]) rank_mapping1[j]=i1++;
                if (i2<two._var.size() && res._var[j]==two._var[i2]) rank_mapping2[j]=i2++;
            }
        }
        //cout<<"rank_mapping1="<<rank_mapping1<<endl;
        //cout<<"rank_mapping2="<<rank_mapping2<<endl;
        
        // main loop
        one.reset();
        two.reset();
        res.reset();
        int rank=0;
        while (rank<res._size) {
            //cout<<"***rank="<<rank<<endl;
            //cout<<"one._current_index="<<one._current_index<<endl;
            //cout<<"two._current_index="<<two._current_index<<endl;
            //cout<<"res._current_index="<<res._current_index<<endl;
            res._data[res._current_index]=one.get()*two.get();
            rank=0;
            one.increment_index(rank_mapping1[rank]);
            two.increment_index(rank_mapping2[rank]);
            while (res.increment_index(rank)==0) {
                rank++;
                if (rank<res._size) {
                    one.increment_index(rank_mapping1[rank]);
                    two.increment_index(rank_mapping2[rank]);
                }
            }
        }
        // return result
        
        return(res);
    }
    
    // affect val to all potential configurations compatible with conf
    void affect(map<string,string> conf,valarray<double> val) {
        if (_size==0) {
            _sum=val;
        }
        // build target conf
        vector<int> target_conf;
        target_conf.resize(_size,-1);
        for (int i=0; i<_size; i++) {
            map<string,string>::iterator it=conf.find(_var[i].name());
            if (it!=conf.end())
                for (int j=0; j<_dim[i]; j++)
                    if (it->second==_var[i]._val[j]) target_conf[i]=j;
        }
        //cout<<"target_conf="<<target_conf<<endl;
        // main loop
        reset();
        int rank=0;
        while (rank<_size) {
            //cout<<"***rank="<<rank<<endl;
            //cout<<"_current_conf="<<_current_conf<<endl;
            bool compatible=true;
            for (int i=0; i<_size; i++) if ((_current_conf[i]!=target_conf[i])&(target_conf[i]>-1)) {
                compatible=false;
                break;
            }
            if (compatible==true)
                _data[_current_index]=val;
            // increment
            rank=0;
            while (increment_index(rank)==0) rank++;
        }
    }
    void affect(map<string,string> conf,double val) {
        valarray<double> vval;
        vval.resize(SIZE,val);
        affect(conf,vval);
    }
    void affect(valarray<double> val){
        map<string,string> conf;
        affect(conf,val);
    }
    void affect(double val){
        valarray<double> vval;
        vval.resize(SIZE,val);
        affect(vval);
    }
    
    // change variable name in a potential
    void rename(map<string,string> sub) {
        map<string,string>::iterator it;
        for (int i=0; i<_size; i++) {
            it=sub.find(_var[i].name());
            if (it!=sub.end())
                _var[i].rename(it->second);
        }
        sort(_var.begin(),_var.end());
    }

};

class bn {
  
    // main members
    vector<variable> _var;
    vector<potential> _K;
    // boolean status
    bool _valid;
    bool _jt;
    bool _bp;
    // redundant information
    int _size;
    map<string,int> _id;
    vector<string> _sid;
    vector<set<int> > _par;
    vector<int> _dim;
    // junction tree members
    vector<set<int> > _adj;
    vector<int> _fillin;
    vector<bool> _eliminated;
    vector<set<int> > _cl;
    vector<set<int> > _sep;
    vector<int> _link;
    vector<int> _cliqueof; // clique membership for all variables
    // belief propagation
    vector<set<int> > _neighbors; // clique neighboring
    vector<potential> _C; // potential of all cliques
    // workspace for computations
    vector<potential> _inward;
    vector<potential> _outward;
    vector<potential> _marginalcl;
    vector<potential> _marginalsep;
    valarray<double> _pevidence;

    
    // set union and intersection
    set<int> union_set(set<int> &s1,set<int> &s2) {
        set<int> res;
        set_union(s1.begin(),s1.end(),s2.begin(),s2.end(),inserter(res, res.begin()));
        return(res);
    }
    set<int> inter_set(set<int> &s1,set<int> &s2) {
        set<int> res;
        set_intersection(s1.begin(),s1.end(),s2.begin(),s2.end(),inserter(res, res.begin()));
        return(res);
    }
    // compute the number of fillin needed to turn the neighboring of i into a clique
    void compute_fillin(int i) {
        _fillin[i]=0;
        if (_adj[i].size()<=1) {
            // clique size 2 or less
            // no fill-in !
        } else {
            // clique of size 3 or more
            set<int> clique=_adj[i];
            clique.insert(i);
            for (set<int>::iterator it=_adj[i].begin();it!=_adj[i].end(); it++) {
                _fillin[i]+=_adj[i].size()-inter_set(clique,_adj[*it]).size();
            }
        }
        _fillin[i]/=2;
    }
    void initialize_fillin() {
        // compute number of fill-ins
        _fillin.resize(_sid.size());
        for (int i=0; i<_sid.size(); i++)
            compute_fillin(i);
    }

    // recursive JT building functions
    int update(set<int> &elim) {
        // loop on elim set
        for (set<int>::iterator it=elim.begin(); it!=elim.end(); it++) {
            // remove *it from adj
            for (set<int>::iterator jt=_adj[*it].begin(); jt!=_adj[*it].end(); jt++) {
                _adj[*jt].erase(_adj[*jt].find(*it));
            }
            // update fillin
            for (set<int>::iterator jt=_adj[*it].begin(); jt!=_adj[*it].end(); jt++) {
                compute_fillin(*jt);
            }
            // update eliminated
            _eliminated[*it]=true;
        }
        // return next elimination variable
        int var=-1;
        int best=INT_MAX;
        for (int i=0; i<_eliminated.size(); i++) {
            if (!_eliminated[i]) {
                if (_fillin[i]<best) {
                    best=_fillin[i];
                    var=i;
                }
            }
            if (best==0) break;
        }
        return(var);
    }
    set<int> remove(int i) {
        // get clique
        set<int> clique=_adj[i];
        clique.insert(i);
        // add fill-in if necessary
        set<int> to_update;
        if (_fillin[i]>0) {
            for (set<int>::iterator it=_adj[i].begin();it!=_adj[i].end(); it++) {
                // fillin update list
                to_update.insert(*it);
                for (set<int>::iterator jt=_adj[*it].begin(); jt!=_adj[*it].end(); jt++) {
                    to_update.insert(*jt);
                }
                // add fill-in
                for (set<int>::iterator jt=clique.begin(); jt!=clique.end(); jt++) {
                    if (*jt!=*it) {
                        _adj[*it].insert(*jt);
                        _adj[*jt].insert(*it);
                    }
                }
            }
        }
        // update fillin
        for (set<int>::iterator it=to_update.begin();it!=to_update.end(); it++)
            compute_fillin(*it);
        // separator set
        set<int> separator;
        set<int> elim;
        elim.insert(i);
        for (set<int>::iterator it=_adj[i].begin();it!=_adj[i].end(); it++) {
            if (inter_set(clique,_adj[*it]).size()!=_adj[*it].size())
                separator.insert(*it);
            else
                elim.insert(*it);
        }
        // add clique and separator
        _cl.push_back(clique);
        _sep.push_back(separator);
        return(elim);
    }

    
public:
    
    // empty constructor
    bn() {
        _size=0;
        _valid=false;
        _jt=false;
        _bp=false;
    };
    
    // add a couple variable potential to the BN
    // check in the process that variable V is included in K
    void add(const variable &V,const potential &K) {
        bool valid=false;
        for (vector<variable>::const_iterator it=K._var.begin(); it!=K._var.end(); it++)
            if (it->_id==V._id) {
                valid=true;
                break;
            }
        if (!valid) {
            cerr<<"error in bn::add(): couple (variable,potential) is not valid"<<endl;
            return;
        }
        // todo
        _var.push_back(V);
        _K.push_back(K);
        _size++;
    };
    
    bool initialize(bool verbose=true) {
        // check for duplicated variables
        for (int i=0; i<_size; i++)
            _id[_var[i]._id]=i;
        if (_id.size()<_size) {
            cout<<"error in bn:initialize(): duplicated variables"<<endl;
            return false;
        }
        // fill redundant (but useful) information: _par,_dim,_sid;
        _par.resize(_size);
        _dim.resize(_size);
        _sid.resize(_size);
        for (int i=0; i<_size; i++) {
            _dim[i]=_var[i].dim();
            _sid[i]=_var[i]._name;
            for (vector<variable>::iterator it=_K[i]._var.begin(); it!=_K[i]._var.end(); it++) {
                // check if the variable is in _id
                map<string,int>::iterator iit=_id.find(it->_id);
                if (iit==_id.end()) {
                    cout<<"error in bn:initialize(): forbidden variable in a potential"<<endl;
                    return false;
                }
                if (iit->second!=i) _par[i].insert(iit->second);
            }
        }
        _valid=true;
        if (verbose) {
            cout<<endl<<"Bayesian network successfully built:"<<endl;
            cout<<"\tnumber of variables = "<<_id.size()<<endl;
            {
                int nedges=0;
                for (vector<set<int> >::iterator it=_par.begin(); it!=_par.end(); it++)
                    nedges+=it->size();
                cout<<"\tnumber of edges = "<<nedges<<endl;
            }
            {
                double lognconfs=0.0;
                for (vector<int>::iterator it=_dim.begin(); it!=_dim.end(); it++)
                    lognconfs+=log((double)*it);
                cout<<"\ttotal complexity ~ 10^"<<(int)(lognconfs/log(10))<<endl<<endl;
            }
        }
        return _valid;
    }
    
    void bndot(string dotfile="bn.dot") {
        if (!_valid) {
            cerr<<"error in bn::bndot(): cannot produce a dotfile from an invalid bn"<<endl;
            return;
        }
        ofstream out(dotfile);
        out<<"/** dot "<<dotfile<<" -Tpdf > bn.pdf && open bn.pdf **/"<<endl;
        out<<"digraph \"BN generated by bpbat\" {"<<endl;
        out<<"nodesep=0.3;"<<endl;
        out<<"center=true;"<<endl;
        out<<"rankdir=TB;"<<endl;
        out<<"node [peripheries=1,style=\"filled\",bgcolor=\"grey\"];"<<endl;
        //out<<"edge [dir=\"none\"];"<<endl;
        for (int i=0; i<_size; i++) {
            out<<_sid[i]<<";"<<endl;
            for (set<int>::iterator it=_par[i].begin(); it!=_par[i].end(); it++) {
                out<<_sid[*it]<<" -> "<<_sid[i];//<<" [label=\" ";
                out<<";"<<endl;
            }
        }
        out<<"}"<<endl;
    }
    
    void jt(bool verbose=true) {
        // allocate workspace
        _fillin.resize(_size,0);
        _eliminated.resize(_size,false);
        _adj.resize(_size);
        // build the moral graph
        for (int id=0; id<_size; id++) {
            //cout<<"id="<<id<<endl;
            for (set<int>::iterator it=_par[id].begin(); it!=_par[id].end(); it++) {
                //cout<<*it<<endl;
                _adj[id].insert(*it);
                _adj[*it].insert(id);
                for (set<int>::iterator jt=_par[id].begin(); jt!=_par[id].end(); jt++) {
                    if (*it!=*jt) {
                        _adj[*it].insert(*jt);
                        _adj[*jt].insert(*it);
                    }
                }
            }
        }
        initialize_fillin();
        // perform variable elimination
        set<int> elim;
        int next=update(elim);
        while (next>-1) {
            //print_adj();
            //cout<<"remove = "<<_sid[next]<<endl;
            // remove elim from graph
            elim=remove(next);
            //cout<<"elim = "; print_set(elim); cout<<endl;
            // display clique and separator
            //cout<<"cl = "; print_set(_cl[_cl.size()-1]); cout<<endl;
            //cout<<"sep = "; print_set(_sep[_sep.size()-1]); cout<<endl;
            // next elimination
            next=update(elim);
            //cout<<"next = "<<_sid[next]<<endl;
        }
        // connect cliques
        _link.resize(_cl.size(),-1);
        for (int i=0; i<_cl.size(); i++) {
            if ((i+1)<_cl.size()) {
                for (int j=i+1; j<_cl.size(); j++) {
                    if (inter_set(_sep[i],_cl[j]).size()==_sep[i].size()) {
                        _link[i]=j;
                        break;
                    }
                }
            }
        }
        // affect variables to cliques
        _cliqueof.resize(_size);
        for (int i=0; i<_size; i++) {
            int j;
            set<int> fam=_par[i];
            fam.insert(i);
            for (j=0; j<_size; j++) {
                if (inter_set(fam,_cl[j]).size()==fam.size())
                    break;
            }
            //cout<<"_cliqueof["<<i<<"]="<<j<<endl;
            _cliqueof[i]=j;
        }
        // report _jt built
        _jt=true;


        if (verbose) {
            cout<<endl<<"Junction tree successively built:"<<endl;
            int complexity=0;
            int tree_width=0;
            for (vector<set<int> >::iterator it=_cl.begin(); it!=_cl.end(); it++) {
                int conf=1;
                if (it->size()>tree_width)
                    tree_width=it->size();
                for (set<int>::iterator jt=it->begin(); jt!=it->end(); jt++)
                    conf=conf*_dim[*jt];
                complexity+=conf;
            }
            cout<<"\tnumber of cliques = "<<_cl.size()<<endl;
            cout<<"\tmax clique size = "<<tree_width<<endl;
            cout<<"\ttotal complexity = "<<complexity<<endl<<endl;
        }
    };
    
    void jtdot(string dotfile="jt.dot") {
        if (!_jt) {
            cerr<<"error in bn::jtdot(): cannot produce a dotfile for a non-built jt"<<endl;
            return;
        }
        ofstream out(dotfile);
        out<<"/** dot "<<dotfile<<" -Tpdf > jt.pdf && open jt.pdf **/"<<endl;
        out<<"digraph \"JT generated by bpbat\" {"<<endl;
        out<<"nodesep=0.3;"<<endl;
        out<<"center=true;"<<endl;
        out<<"rankdir=BT;"<<endl;
        out<<"node [peripheries=1,style=\"filled\",bgcolor=\"grey\"];"<<endl;
        out<<"edge [dir=\"none\"];"<<endl;
        // all cliques
        for (int i=0; i<_cl.size(); i++) {
            out<<"C"<<i<<" [label=\"C"<<(i+1)<<": ";
            for (set<int>::iterator it=_cl[i].begin(); it!=_cl[i].end(); it++) {
                out<<_sid[*it];
                if (_cliqueof[*it]==i) out<<"*";
                out<<" ";
            }
            out<<"\"];"<<endl;
            if (_link[i]>-1) {
                out<<"C"<<i<<" -> C"<<_link[i]<<" [label=\" ";
                for (set<int>::iterator it=_sep[i].begin(); it!=_sep[i].end(); it++)
                    out<<_sid[*it]<<" ";
                out<<"\"];"<<endl;
            }
        }
        out<<"}"<<endl;
    }

    void bp(bool verbose=true) {
        if (!_jt) {
            cerr<<"error in bn::bp(): cannot perform bp without jt"<<endl;
            return;
        }
        // fills neighboring in the jt
        _neighbors.resize(_cl.size());
        for (int i=0; i<_cl.size(); i++) {
            if (_link[i]>-1) {
                _neighbors[i].insert(_link[i]);
                _neighbors[_link[i]].insert(i);
            }
        }
        // fills _C from _K and _cliqueof
        potential neutral;
        neutral.affect(1.0);
        _C.resize(_cl.size(),neutral);
        //cout<<"ok before _C"<<endl;
        for (int i=0; i<_size; i++) {
            //cout<<"***i="<<i<<endl;
            //cout<<"_C["<<_cliqueof[i]<<"]:"<<endl<<_C[_cliqueof[i]]<<endl;
            //cout<<"_K["<<i<<"]:"<<endl<<_K[i]<<endl;
            _C[_cliqueof[i]]=_C[_cliqueof[i]]*_K[i];
        }
        //cout<<"ok before inward"<<endl;
        for (int i=0; i<_cl.size(); i++)
            //cout<<"_C["<<(i+1)<<"]:"<<endl<<_C[i]<<endl;
        // inward
        _inward.resize(_cl.size());
        for (int i=0; i<_cl.size(); i++) {
            potential res=_C[i];
            for (set<int>::iterator it=_neighbors[i].begin(); it!=_neighbors[i].end(); it++)
                if (*it<i) res=res*_inward[*it];
            set<variable> sep;
            for (set<int>::iterator it=_sep[i].begin(); it!=_sep[i].end(); it++)
                sep.insert(_var[*it]);
            _inward[i]=res.marginal(sep);
            _inward[i].rescale();
            //cout<<"_inward["<<i<<"]:"<<endl<<_inward[i]<<endl;
        }
        // evidence already available
        _pevidence=_inward[_cl.size()-1]._logscale;
        _pevidence+=log(_inward[_cl.size()-1].sum());
        //cout<<fixed<<"log(pev)="<<_pevidence<<endl;
        // outward
        _outward.resize(_cl.size(),neutral);
        if (_cl.size()>1)
            _outward[_cl.size()-1]=neutral;
            for (int i=_cl.size()-2; i>=0; i--) {
                potential res;
                int j=_link[i];
                if (j>-1) {
                    res=_C[j]*_outward[j];
                    for (set<int>::iterator it=_neighbors[j].begin(); it!=_neighbors[j].end(); it++) {
                        if (*it<j && *it!=i) res=res*_inward[*it];
                    }
                }
                set<variable> sep;
                for (set<int>::iterator it=_sep[i].begin(); it!=_sep[i].end(); it++)
                    sep.insert(_var[*it]);
                _outward[i]=res.marginal(sep);
                _outward[i].rescale();
            }
        //cout<<"outward OK"<<endl;
        
        // marginals and checks
        //_marginalsep.resize(_cl.size());
        //for (int i=0; i<_cl.size(); i++) {
        //    _marginalsep[i]=_inward[i]*_outward[i];
        //}
        _marginalcl.resize(_cl.size());
        for (int i=0; i<_cl.size(); i++) {
            _marginalcl[i]=_C[i]*_outward[i];
            for (set<int>::iterator it=_neighbors[i].begin(); it!=_neighbors[i].end(); it++) {
                if (*it<i) _marginalcl[i]=_marginalcl[i]*_inward[*it];
            }
            //cout<<marginalcl[i].sum()<<endl;
            //double relerr=abs((_marginalcl[i].sum()-_pevidence)/_pevidence).max();
            //if (abs(_pevidence).max()<1e-10)
            //    relerr=abs(_marginalcl[i].sum()-_pevidence).max();
            //cout<<"relerr="<<relerr<<endl;
            //if (relerr>1e-10)
            //    cerr<<"error in marginalcl["<<i<<"]: relerr on pevidence = "<<relerr<<endl;
        }
        

        _bp=true;
        if (verbose) {
            cout<<endl<<"Belief propagation successful:"<<endl;
            cout<<"\tlog(evidence probability) = "<<fixed<<_pevidence<<endl<<endl;
        }

    }
    
    potential posterior(variable V) {
        potential res;
        if (!_bp) {
            cerr<<"error in bn::posterior(): cannot compute posterior without bp"<<endl;
            return res;
        }
        int i;
        for (i=0; i<_size; i++)
            if (V._id==_var[i]._id) break;
        if (i>=_size) {
            cerr<<"error in bn::posterior(): unknown variable provided."<<endl;
            return res;
        }
        res=_marginalcl[_cliqueof[i]].marginal({V});
        res._logscale-=_pevidence;
        //potential ev;
        //ev.affect(-_pevidence);
        //res=res*ev;
        //res.logscale(false);
        return res;
    }
        
    valarray<double> pevidence() {
        valarray<double> res(0.,SIZE);
        if (!_bp) {
            cerr<<"error in bn::pevidence(): cannot compute pevidence without bp"<<endl;
            return(res);
        }
        res=_pevidence;
        return(res);
    }

};

#endif

