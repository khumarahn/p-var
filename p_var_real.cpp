// p-variation calculus for piecewise monotone functions
// Author: Vygantas Butkus <Vygantas.Butkus@gmail.com>
// Please do not hesitate to contact me in any question.

#include "p_var_real.h"

namespace p_var_real {

	// -------------------------------- definitions of types  ---------------------------------- //

	// p-variation point. An object with necessary info.
	struct pvpoint {
		int id;
		double val;
		double pvdiff;
	};

	typedef std::list<pvpoint> PrtList;    // the list of p-variation point
	typedef PrtList::iterator it_PrtList;  // the iterator of element in list of p-variation point

	// p-variation temporary points, used in specific calculations -
	// it saves an extra value
	struct pvtemppoint{
		it_PrtList it;
		double ev;
	};

	// last iterator of the list. It is iterator to `obj.end()-1`
	template <class T>
		typename T::iterator last(T& obj){
			typename T::iterator it = obj.end();
			--it;
			return(it);
		}

	// the difference used in p-variation, i.e. the abs power of diff.
	double pvar_diff(double diff, double p){
		return std::pow(std::abs(diff), p);
	}

	// finds change point of the vector and put it in list
	PrtList ChangePoints(const NumericVector& x){

		// Main principle:
		// if point pt[i] is in increasing interval and pt[i]>pt[i+1], then
		// pt[i] is change point (and vise versa).

		int dir = 0;
		int n = x.size();

		pvpoint pvp;
		pvp.id = 0;
		PrtList out (1, pvp); // the first and last points are always included by definition.

		for(int i = 1; i < n; ++i) {
			if(x[i-1]<x[i]){
				if(dir<0){
					pvp.id = i-1;
					out.push_back (pvp);
				}
				dir = 1;
			}
			if(x[i-1]>x[i]){
				if(dir>0){
					pvp.id = i-1;
					out.push_back (pvp);
				}
				dir = -1;
			}
		}

		pvp.id =n-1;
		out.push_back (pvp); // the first and last points are always included by definition.

		return(out);
	}

	// finds(updates) all necessary attributes of `prt` list. `prt` must have good ids.
	void prepare_prt(const NumericVector& x, std::list<pvpoint>& prt,  const double& p){

		PrtList::iterator it1_prt, it2_prt;
		it1_prt = it2_prt = prt.begin();
		++it2_prt;

		// getting first element
		(*it1_prt).val = x[(*it1_prt).id];
		(*it1_prt).pvdiff = 0;

		// getting all other elements
		for ( ; it2_prt != prt.end(); it1_prt++, it2_prt++){
			(*it2_prt).val = x[(*it2_prt).id];
			(*it2_prt).pvdiff = pvar_diff(x[(*it1_prt).id] - x[(*it2_prt).id], p);
		}
	}

	// sequentially checks small intervals of length d
	void CheckSmallIntervalsOnce(PrtList& prt, const double& p,  const int& d){

		// Main principle:
		// if |pt[i] - pt[i+ d]|^p > sum_{j={i+1}}^d   |pt[j] - pt[j-1]|^p
		// then all middle points (i.e. p[j], j=i+1,...,i+d-1) are insignificant

		int dcount = 0;
		double csum = 0;
		double fjoinval;

		PrtList::iterator it1_prt, it2_prt;
		it1_prt = it2_prt = prt.begin();
		++it2_prt;

		for ( ; it2_prt != prt.end(); it2_prt++){
			++dcount;
			csum += (*it2_prt).pvdiff ;
			if(dcount==d){
				fjoinval = pvar_diff((*it1_prt).val - (*it2_prt).val, p);
				++it1_prt;  // in this stage it1_prt is always significant
				if(csum < fjoinval){ // mid points are insignificant, delete all
					dcount = 0;
					csum = 0;
					it1_prt = prt.erase(it1_prt, it2_prt);
					(*it2_prt).pvdiff = fjoinval;
				}else{ // all significant, move one step foward forward
					csum -= (*it1_prt).pvdiff ;
					--dcount;
				}
			}
		}
	}

	// checks small intervals of up till length dn.
	// After this function, all points are significant in any small interval (i.e. interval with length no greater dn).
	void CheckSmallIntervals(PrtList& prt, const double& p,  const int& dn){

		// Main principle:
		// apply CheckSmallIntervalsOnce starting form d=3 (because 3 is the minimal length worth checking)
		// If there was no change, apply CheckSmallIntervalsOnce with d=d+2 (because insignificant points goes only in pears)
		// If there was change start from d=3 again (because `prt` changed, therefore, we are not sure if all smaller intervals are good).

		int LastSize = 0;
		int CurSize = prt.size();
		int d = 3;

		while((LastSize!=CurSize) & (CurSize>3) & (d<=dn)){
			d = 3;
			LastSize = CurSize;
			CheckSmallIntervalsOnce(prt, p, d);
			CurSize = prt.size();
			while((LastSize==CurSize) & (CurSize>d+2) & (d<dn)){
				d = d + 2;
				LastSize = CurSize;
				CheckSmallIntervalsOnce(prt, p, d);
				CurSize = prt.size();
			}
		}
	}


	// merge two intervals ([a, v] and [v, b]) which are known to be good.
	void Merge2GoodInt(PrtList& prt,  const double& p, it_PrtList a, it_PrtList v, it_PrtList b){

		// Main principle:
		// 1. Find potential points in intervals [a,v) and (v, b]
		//    (i.e. the points that could make a new f-joint with any point form opposite interval).
		//    Those points are find using cummin and cummac starting from v.
		//     Some points might be dropped out before actual checking, but experiment showed, that it is not worthwhile.
		// 2. Sequentially check all possible joints. If any increase is detected, then all middle points are insignificant.

		if (a==v or v==b) return ; // nothing to calculate, exit the procedure.

		double amin, amax, bmin, bmax, ev, balance, maxbalance, jfoin, takefjoin;
		it_PrtList prt_it, prt_ait, prt_bit;
		std::list<pvtemppoint> av, vb;
		std::list<pvtemppoint>::iterator ait, bit, tit, tait, tbit, bitstart;
		pvtemppoint pvtp;

		// 1. ### Find potential points

		// --- in interval [a,v) (starting from v).
		ev = 0;
		prt_it = v;
		amin = amax= (*v).val;
		while(prt_it!=a){
			ev += (*prt_it).pvdiff;
			--prt_it;
			if((*prt_it).val>amax){
				amax=(*prt_it).val;
				pvtp.it = prt_it;
				pvtp.ev = ev;
				av.push_back (pvtp);
			}
			if((*prt_it).val<amin){
				amin=(*prt_it).val;
				pvtp.it = prt_it;
				pvtp.ev = ev;
				av.push_back (pvtp);
			}
		}
		// printList(av, "av :");

		// --- in interval (v,b] (starting from v).
		ev = 0;
		prt_it = v;
		bmin = bmax = (*v).val;
		while(prt_it!=b){
			++prt_it;
			ev += (*prt_it).pvdiff;
			if((*prt_it).val>bmax){
				bmax=(*prt_it).val;
				pvtp.it = prt_it;
				pvtp.ev = ev;
				vb.push_back (pvtp);
			}
			if((*prt_it).val<bmin){
				bmin=(*prt_it).val;
				pvtp.it = prt_it;
				pvtp.ev = ev;
				vb.push_back (pvtp);
			}
		}
		// printList(vb, "vb :");

		// 2. ### Sequentially check all possible joints: finding the best i,j \in [a, v)x(v,b] that could be joined
		takefjoin = 0;
		maxbalance = 0;
		for(ait=av.begin(); ait!=av.end(); ait++){
			for(bit=vb.begin(); bit!=vb.end(); bit++){
				// std::cout <<  (*(*ait).it).id << " - " << (*(*bit).it).id << ":\n";
				jfoin = pvar_diff( (*(*ait).it).val - (*(*bit).it).val, p );
				balance = jfoin - (*bit).ev - (*ait).ev ;
				if (balance>maxbalance){
					maxbalance = balance;
					takefjoin = jfoin;
					tait = ait;
					tbit = bit;
				}
			}
		}

		// if we found any point, join it by erasing all middle points
		if(maxbalance>0){
			// joining:
			prt_ait = (*tait).it;
			++prt_ait;
			prt_it = prt.erase(prt_ait, (*tbit).it);
			(*prt_it).pvdiff = takefjoin;
		}
	}

	// Modifies prt to become the partition of p-variation by merging all small in [a, b]
	void PvarByMerging(PrtList& prt,  const double& p, PrtList::iterator a, PrtList::iterator b, int LSI=2){

		// Main principle:
		// 1. find all the intervals that should be merged.
		// 2. Apply merging by pears of interval. Repeat it until all intervals are merged.

		it_PrtList it1_prt, it2_prt, v, it, it2;
		std::list<it_PrtList> IterList;
		std::list<it_PrtList>::iterator it_IterList, a_IL, v_IL, b_IL;

		// 1. ### Finding all the intervals that will be merged
		it = a;
		int count = 0;
		while(it!=b){
			if(count % LSI == 0){
				IterList.push_back (it);
			}
			++count;
			++it;
		}
		IterList.push_back (it);

		// ### 2. Apply merging by pears of interval until everything is merged.
		while(IterList.size()>2){
			a_IL = v_IL = b_IL = IterList.begin();
			++v_IL;
			++b_IL;
			++b_IL;
			while((b_IL!=IterList.end()) & (v_IL!=IterList.end())){
				Merge2GoodInt(prt, p, *a_IL, *v_IL, *b_IL);
				a_IL = IterList.erase(v_IL);
				v_IL = b_IL = a_IL;
				++v_IL;
				++b_IL;
				++b_IL;
			}
		}
	}

	// p-variation calculation (in C++)
	double pvar(const NumericVector& x, double p) {

		const int LSI=3;

		// the length of x
		if (x.size() < 1) {
			return 0;
		}
		if (x.size() == 2) {
			return std::pow(std::abs(x[0]-x[1]), p);
		}

		// ##### Program itself:

		PrtList::iterator it_prt;
		PrtList prt = ChangePoints(x);
		prepare_prt(x, prt, p);

		CheckSmallIntervals(prt, p, LSI);
		PvarByMerging(prt, p, prt.begin(), last(prt), LSI+1);

		// output:
		double pvalue=0;;
		NumericVector partition(prt.size());
		int i = 0;
		for(it_prt=prt.begin(); it_prt!=prt.end(); it_prt++){
			pvalue += (*it_prt).pvdiff;
			partition[i] = (*it_prt).id + 1;
			++i;
		}

		// pvalue = std::pow(pvalue, 1/p);

		return pvalue;
	}

} // namespace
