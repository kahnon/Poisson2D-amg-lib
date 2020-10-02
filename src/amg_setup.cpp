#include <vector>
#include <tuple>
#include <set>
#include <algorithm>
#include <iostream>
#include <math.h>

#include "amg_setup.h"
#include "matrix.h"

#define USE_SECOND_PASS 0

void coarsen_matrix(Matrix& matrix, double thresh, double beta){
  //std::set combining weight and index of each grid point which is used for selection of coarse points
  //use std::set which is ordered with passed ordering function
  typedef std::set<std::tuple<int,double>,bool(*)(const std::tuple<int,double>& t1, const std::tuple<int,double>& t2)> WeightSet;

  typedef WeightSet::iterator iterator_tp;

  WeightSet weights([](const std::tuple<int,double>& t1, const std::tuple<int,double>& t2){
      double val1 = std::get<double>(t1);
      double val2 = std::get<double>(t2);
      int ind1 = std::get<int>(t1);
      int ind2 = std::get<int>(t2);
      //otherwise sort by descending weight primarily and ascending index secondarily
      if(val1!=val2) return val1 > val2;
      else return ind1 < ind2;
    }
  );


  //array associating iterator and array_index of each point for fast element access
  std::vector<iterator_tp> weights_it(matrix.size());
  
  //necessary functions for coarsening algorithm
  auto calc_weight = [&](int i){
    double wt = 0;
    auto& ln = matrix[i];

    if(ln.status==UNDECIDED){
      for(auto& cn : ln.conns){
	int ind = std::get<int>(cn);
	if(matrix[ind].status==SFINE && std::get<ConnType>(cn)==STRONG) wt += 2;
	else if(matrix[ind].status==UNDECIDED && std::get<ConnType>(cn)==STRONG) wt += 1;
      }
    }
    return wt;
  };

  //returns coupling strength of j to Coarse_points in C(i) and from i to j
#if USE_SECOND_PASS
  auto calc_coupling_strength = [&] (int j, int i){
    //coupling fro mj to C(i)
    double strength_j_to_ci = 0;

    //not sure if diagonal values for i and j should be considered
    double max_j_entry = fabs(matrix[j].val);
    double max_i_entry = fabs(matrix[i].val);

    for(auto& cn : matrix[j].conns){
      double val = std::get<double>(cn);
      if(fabs(val) > max_j_entry) max_j_entry = fabs(val);
    }

    //add entries in coarse neighborhood of 
    for(auto& cn : matrix[i].conns){
      int ind = std::get<int>(cn);
      double val = std::get<double>(cn);
      //find max i entry
      if(fabs(val) > max_i_entry) max_i_entry = fabs(val);

      //calc coupling strength
      if(std::get<ConnType>(cn) == STRONG && matrix[ind].status==COARSE){
	strength_j_to_ci -= find_matrix_val(matrix,j,ind);
      }
    }
    assert(max_j_entry != 0 && max_i_entry != 0 && "max entry is 0 in calc_coupling_strength");
    strength_j_to_ci /= max_j_entry;

    //now calc strength of coupling from i to j
    double strength_i_to_j = find_matrix_val(matrix,i,j)/max_i_entry;

    return std::make_pair(strength_j_to_ci,strength_i_to_j);
  };
#endif


  //remove points in list to_update
  auto erase_weights = [&](std::vector<int> to_erase){
    for(auto ind : to_erase){ 
      iterator_tp& it = weights_it[ind];
      weights.erase(it);
      it = weights.end();
    }
  };

  //update points in list to_update
  auto update_weights = [&](std::vector<int> to_update){
    //erase old tuple, update weight of index an insert new tuple into weights, then update iterator
    for(auto ind : to_update){
      iterator_tp& it = weights_it[ind];
      double wt = calc_weight(ind);

      weights.erase(it);
      it = weights.insert(std::make_tuple(ind,wt)).first;
    }
  };

  auto add_coarse_point = [&](int i){
    //std::cout<<"adding point "<<std::get<int>(*weights.begin())<<"  with weight "<<std::get<double>(*weights.begin())<<"  remaining before erasing current "<<weights.size()<<std::endl;
    std::vector<int> to_update;
    std::vector<int> to_erase;
    auto& coarse_point = matrix[i];

    //for second pass coarsening
    if(coarse_point.status != UNDECIDED){
      coarse_point.status=COARSE;
    //first pass coarsening
    }else{
      coarse_point.status=COARSE;
      to_erase.push_back(i);

      //first mark all strongly connected points as SFINE, then add them to erase them from vector
      for(auto& coarse_conn : coarse_point.conns){
	int ind = std::get<int>(coarse_conn);
	auto& fine_point = matrix[ind];
	if(fine_point.status == UNDECIDED && std::get<ConnType>(coarse_conn) == STRONG){
	  fine_point.status = SFINE;
	  if(weights_it[ind] != weights.end()) to_erase.push_back(ind);
	}
      }
      //erase those points from std::set
      erase_weights(to_erase);
      //update the weights of weakly connected neighbors
      //i's conncetions as to_update if they are undecided. 
      for(auto& coarse_conn : coarse_point.conns){
	int f_ind = std::get<int>(coarse_conn);
	auto& fine_point = matrix[f_ind];
	if(weights_it[f_ind] != weights.end()) to_update.push_back(f_ind);

	for(auto& fine_conn : fine_point.conns){
	  int ind = std::get<int>(fine_conn);
	  if(matrix[ind].status == UNDECIDED){ 
	    if(weights_it[ind] != weights.end()){
	      if(std::find(to_update.begin(),to_update.end(),ind) == to_update.end()) to_update.push_back(ind);
	    }
	  }
	}
      }
      //output for debugging
#if 0
      for(auto& ind : to_update)
	std::cout<<"   updating point "<<ind<<"  old weight "<<std::get<double>(*weights_it[ind])<<"  new_weight "<<calc_weight(ind)<<std::endl;
      for(auto& ind : to_erase)
	std::cout<<"   erasing  point "<<ind<<std::endl;
#endif
      update_weights(to_update);
    }
  };

  //**********ACTUAL COARSENING ALGORITHM*************
  //first set transpose connections
  set_matrix_transpose_connections(matrix,thresh);

  //calculate weights and insert them into the set which takes care of the order
  for(auto i=0u;i<weights_it.size();++i){
    double wt = calc_weight(i);
    std::pair<iterator_tp,bool> it_pair = weights.insert(std::make_tuple(i,wt));
    weights_it[i]=it_pair.first;
    if(it_pair.second == false){ 
      std::cout<<"point "<<i<<" with weight "<<wt<<" not inserted into set!"<<std::endl;
      weights_it[i]=weights.end();
    }
  }

  //erase 0 weight entries and mark them as boundary entries
  for(auto& it:weights_it){
    if(std::get<double>(*it)==0){
      matrix[std::get<int>(*it)].status = BOUNDARY;
      weights.erase(it);
      it = weights.end();
    }
  }

  //first pass of coarsening
  //creates maximal set in the sense that no coarse point strongly depends on other coarse point
  do{
    add_coarse_point(std::get<int>(*weights.begin()));
  }while(!weights.empty());

#if USE_SECOND_PASS
  //second pass of coarsening ensures that strong F-F connections require a common coarse point
  for(auto i=0;i<matrix.size();++i){
    auto& ln = matrix[i];
    if(ln.status==UNDECIDED){ 
      std::cout<<"still undecided in second pass of coarsening. This should not happen!"<<std::endl;
      exit(-1);
    }
    else if(ln.status==COARSE) continue;
    //only add additional coarse points near fine points
    else{
      int coarse_point_added=0;

      for(auto& cn : ln.conns){
	int ind = std::get<int>(cn);
	if(matrix[ind].status==COARSE || std::get<ConnType>(cn)!=STRONG) continue;

	auto coupling_strength = calc_coupling_strength(ind,i);
	if(coupling_strength.first <= beta* coupling_strength.second){
	  if(coarse_point_added == 0){
	    std::cout<<"added point "<<ind<<" via second pass at i="<<i<<std::endl;
	    matrix[ind].status=COARSE;
	    coarse_point_added=ind;
	  }else{
	    std::cout<<"added point "<<i<<" via second pass"<<std::endl;
	    matrix[i].status=COARSE;
	  }
	}
      }
    }
  }
#endif
}


//creates interpolation from coarser level to the level of matrix
//using direct interpolation as defined by K. Stueben
//direct interpolation is good for zero row sum matrices where each fine point needs at least 1 strong coarse neighbor connection
// means no point is (WFINE)
Matrix create_interpolation_matrix(Matrix& matrix, double thresh){
  Matrix interpolation;
  //mapping of coarse points between old and new grid
  std::vector<int> coarse_map(matrix.size(),-2);
  int coarse_ctr=0;
  for(auto i=0u;i<matrix.size();++i){
    if(matrix[i].status == COARSE){
      coarse_map[i]=coarse_ctr++;
    }
  }
  //calc weight for one interpolation matrix element. i is row index, j is column index (if i!=j -> goes into conns)
  //direct interpolation, is usually rather bad
#if 0
  //uncommented to eliminate compiler warning
  //not sure what I used this for
  auto calc_interpolation_weights_direct = [&](int i){
    std::vector<std::tuple<int,double,ConnType>> weights;
    const line& ln = matrix[i];
    double sum_pos_interpolation_points = 0;
    double sum_neg_interpolation_points = 0;
    double alpha = 0;
    double beta = 0;

    for(auto& cn : ln.conns){
      double val = std::get<double>(cn);
      int ind = std::get<int>(cn);

      if(val>0){
	if(matrix[ind].status == COARSE && std::get<ConnType>(cn)==STRONG) sum_pos_interpolation_points += val;
	beta += val;
      }else if(val<0){
	if(matrix[ind].status == COARSE && std::get<ConnType>(cn)==STRONG) sum_neg_interpolation_points += val;
	alpha += val;
      }
    }

    if(sum_neg_interpolation_points!=0) alpha /= sum_neg_interpolation_points;
    else alpha=0;
    if(sum_pos_interpolation_points!=0) beta /= sum_pos_interpolation_points;
    else beta=0;

    assert((alpha != 0 || beta != 0 ) && "no interpolation points to calc_interpolation_weights");

    for(auto k=0u;k<ln.conns.size();++k){
      auto& cn = ln.conns[k];
      int ind = std::get<int>(cn);
      double val = std::get<double>(cn);

      //the conns are given in the coordinates of the coarser matrix
      if(matrix[ind].status == COARSE && std::get<ConnType>(cn)==STRONG){
	double wt=0;
	if(val<0){
	  wt = -1 * alpha * val / ln.val;
	}else if(val>0){
	  wt = -1 * beta * val / ln.val;
	}
	//do not include weights near zero
	if(fabs(wt)>1e-6) weights.emplace_back(coarse_map[ind],wt,STRONG);
      }
    }
    return line(i,0,SFINE,weights);
  };
#endif

  //calculate coefficient necessary for improved interpolation, formula
  //i=finer grid index, j=coarser grid index
  //c_ij = sum_(k in FINE(i)) ( matrix(i,k)*matrix(k,j) / ((sum_(l in COARSE(i)) (matrix(k,l))) + matrix(k,i))
  auto calc_int_coefficient = [&matrix](int i, int j){
    double coeff = 0;
    for(auto& cn : matrix[i].conns){
      int k = std::get<int>(cn);
      if((matrix[k].status==SFINE || matrix[k].status==WFINE)){
	double k_l_sum = 0;
	for(auto& cn_coarse : matrix[i].conns){
	  int l = std::get<int>(cn_coarse);
	  if(matrix[l].status==COARSE) k_l_sum += find_matrix_val(matrix,k,l);
	  std::cout<<"        k="<<k<<"  l="<<l<<"  sum="<<k_l_sum<<std::endl;
	}
	assert((k_l_sum + find_matrix_val(matrix,k,i)) != 0 && "in create_interpolation_matrix::calc_int_coeff division by 0!"); 
	std::cout<<"     k_l_sum="<<k_l_sum<<"  matrix_val_"<<k<<"_"<<i<<"="<<find_matrix_val(matrix,k,i)<<std::endl;
	coeff += (find_matrix_val(matrix,i,k) * find_matrix_val(matrix,k,j)) / (k_l_sum + find_matrix_val(matrix,k,i));
      }
    }
    return coeff;
  };

  //improved interpolation formula working also for unstructured meshes
  auto calc_interpolation_weights = [&](int i){
    std::cout<<"weights for index "<<i<<":"<<std::endl;
    std::vector<std::tuple<int,double,ConnType>> weights;
    const line& ln = matrix[i];
    assert((ln.status==SFINE || ln.status==WFINE) && "calc'ing interpolation weight for non-FINE point");
    //first calculate coefficient c_ii since it is needed later
    double coeff_ii = calc_int_coefficient(i,i);
    double val_ii   = ln.val;

    for(auto& cn : ln.conns){
      int j = std::get<int>(cn);
      //interpolate from coarse points
      if(matrix[j].status == COARSE){
	double wt = -1*(find_matrix_val(matrix,i,j) + calc_int_coefficient(i,j));
	std::cout<<"  wt="<<wt<<" val_ii="<<val_ii<<" coeff_ii="<<coeff_ii<<std::endl;
	wt /= (val_ii + coeff_ii);
	if(fabs(wt)>1e-6) weights.emplace_back(coarse_map[j],wt,STRONG);
      }
    }
    std::cout<<std::endl;
    return line(i,0,SFINE,weights);
  };

  /******************INTERPOLATION MATRIX CREATION**********************/
  set_matrix_connections(matrix,thresh);

  for(auto i=0u;i<matrix.size();++i){
    const line& ln = matrix[i];
    if(ln.status == UNDECIDED){
      //this really should not happen
      std::cout<<"trying to interpolate a not properly coarsened matrix!"<<std::endl;
      std::cout<<"ind="<<ln.ind<<"  val="<<ln.val<<std::endl<<"neighbors:"<<std::endl;
      for(auto& cn:ln.conns){
	int ind=std::get<int>(cn);
	double val=std::get<double>(cn);
	std::cout<<"   ind="<<ind<<"  val="<<val<<"  connection_type="<<std::get<ConnType>(cn)<<"  stat="<<matrix[ind].status<<std::endl;
      }
      exit(-1);
    }else if(ln.status == COARSE){
      std::vector<std::tuple<int,double,ConnType>> wt;
      wt.emplace_back(coarse_map[i],1,STRONG);
      interpolation.emplace_back(i,0,COARSE,wt);
    }else if(ln.status == SFINE){
      interpolation.emplace_back(calc_interpolation_weights(i));
    }else if(ln.status == WFINE){
      //no point should be WFINE for this type of interpolation
      interpolation.emplace_back(i,0,WFINE);
    }else{
      //should only be boundaries
      interpolation.emplace_back(i,0,BOUNDARY);
    }
  }
  std::cout<<std::endl<<std::endl<<std::endl;
  return interpolation;
}



Matrix create_coarse_matrix(const Matrix& matrix, const Matrix& int_matrix){
  assert(matrix.size() && matrix.size()==int_matrix.size() && "matrix dimensions dont match to create coarse matrix");
  //map of coarse points between old and new grid
  std::vector<int> coarse_map(matrix.size(),-2);
  int coarse_ctr=0;
  for(auto i=0u;i<matrix.size();++i){
    if(matrix[i].status == COARSE){
      coarse_map[i]=coarse_ctr++;
    }
  }

  Matrix coarse_matrix(coarse_ctr);
  for(auto i=0u;i<coarse_matrix.size();++i) coarse_matrix[i].ind=i;

  for(auto k=0u;k<matrix.size();++k){
    const line& k_ln=matrix[k];
    if(k_ln.status==BOUNDARY) continue;

    std::vector<int> k_coarse;
    std::vector<int> l_list(1,k);

    for(auto& cn : k_ln.conns){
      int ind=std::get<int>(cn);
      if(std::get<ConnType>(cn)==STRONG){
	l_list.push_back(ind);
	if(matrix[ind].status==COARSE){  
	  k_coarse.push_back(ind);
	  //std::cout<<"k="<<k<<" added ind="<<ind<<std::endl;
	}
      }
    }

    for(auto& l : l_list){
      std::vector<int> coarse_list(k_coarse);

      //add l if it is Coarse and not yet added
      if((std::find(coarse_list.begin(),coarse_list.end(),l) == coarse_list.end()) && matrix[l].status==COARSE) 
	coarse_list.push_back(l);

      for(auto& cn:matrix[l].conns){
	int ind = std::get<int>(cn);
	if(std::get<ConnType>(cn)==STRONG && matrix[ind].status==COARSE){ 
	  if(std::find(coarse_list.begin(),coarse_list.end(),ind)==coarse_list.end()){
	    coarse_list.push_back(ind);
	    //std::cout<<"   k="<<k<<" l="<<l<<" added ind="<<ind<<std::endl;
	  }
	}
      }

      for(auto i:coarse_list){
	for(auto j:coarse_list){
	  int i_new=coarse_map[i];
	  int j_new=coarse_map[j];
	  double to_add=0;

	  //calc value to add to new matrix element (i_new,j_new)
	  to_add = find_matrix_val(int_matrix,k,i_new)*find_matrix_val(int_matrix,l,j_new)*find_matrix_val(matrix,k,l);
	  //std::cout<<"added "<<to_add<<" to element via indices i="<<i_new<<" j="<<j_new<<" k="<<k<<" l="<<l<<"  "<<find_matrix_val(int_matrix,k,i_new)<<"  "<<find_matrix_val(int_matrix,l,j_new)<<"  "<<find_matrix_val(matrix,k,l)<<std::endl;

	  if(fabs(to_add)>1e-6){
	    //std::cout<<"  added "<<to_add<<" to element via indices i="<<i_new<<" j="<<j_new<<" k="<<k<<" l="<<l<<std::endl;
	    line& new_ln = coarse_matrix[i_new];

	    //if diag_element just add to val, otherwise to the right conn
	    if(i_new==j_new) coarse_matrix[i_new].val+=to_add;
	    else{
	      //if cn exists add there, otherwise create
	      bool added=false;
	      for(auto& cn:new_ln.conns){
		if(std::get<int>(cn) == j_new){
		  std::get<double>(cn) += to_add;
		  added=true;
		}
	      }
	      if(!added){
		new_ln.conns.emplace_back(j_new,to_add,WEAK);
	      }
	    }
	  }
	}
      }
    }
  }

  return coarse_matrix;
}


