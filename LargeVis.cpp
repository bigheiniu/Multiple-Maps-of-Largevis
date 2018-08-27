#include "LargeVis.h"
#include <map>
#include <float.h>
#include<time.h>
#include <numeric>
#include <cmath>
//TODO: use rabbit to support distributed
//problem from pi is greater than 1 => 同步的问题


LargeVis::LargeVis()
{
	vec = vis = weights = prob = NULL;
	knn_vec = old_knn_vec = NULL;
	annoy_index = NULL;
	head = alias = NULL;
	neg_table = NULL;
}

const gsl_rng_type *LargeVis::gsl_T = NULL;
gsl_rng *LargeVis::gsl_r = NULL;

void LargeVis::clean_model()
{
	if (vis) delete[] vis;
	if (prob) delete[] prob;
	if (weights) delete[] weights;
	if (knn_vec) delete[] knn_vec;
	if (old_knn_vec) delete[] old_knn_vec;
	if (annoy_index) delete annoy_index;
	if (neg_table) delete[] neg_table;
	if (alias) delete[] alias;

	weights = vis = prob = NULL;
	knn_vec = old_knn_vec = NULL;
	annoy_index = NULL;
	neg_table = NULL;
	alias = NULL;

	edge_count_actual = 0;
	neg_size = 1e8;
}

void LargeVis::clean_graph()
{
	if (head) { delete[] head; head = NULL; }

	n_edge = 0;
	next.clear(); edge_from.clear(); edge_to.clear(); reverse.clear(); edge_weight.clear(); names.clear();
}

void LargeVis::clean_data()
{
	if (vec) { delete[] vec; vec = NULL; }
	clean_graph();
}

void LargeVis::load_from_file(char *infile)
{
	clean_data();
	FILE *fin = fopen(infile, "rb");
	if (fin == NULL)
	{
		printf("\nFile not found!\n");
		return;
	}
    	printf("Reading input file %s ......", infile); fflush(stdout);
	fscanf(fin, "%lld%lld", &n_vertices, &n_dim);
	vec = new real[n_vertices * n_dim];
	for (long long i = 0; i < n_vertices; ++i)
	{
		for (long long j = 0; j < n_dim; ++j)
		{
			fscanf(fin, "%f", &vec[i * n_dim + j]);
		}
	}
	fclose(fin);
	printf(" Done.\n");
	printf("Total vertices : %lld\tDimension : %lld\n", n_vertices, n_dim);
}

void LargeVis::load_from_data(real *data, long long n_vert, long long n_di)
{
	clean_data();
	vec = data;
	n_vertices = n_vert;
	n_dim = n_di;
	printf("Total vertices : %lld\tDimension : %lld\n", n_vertices, n_dim);
}

void LargeVis::load_from_graph(char *infile)
{
	clean_data();
	char *w1 = new char[1000];
	char *w2 = new char[10000];
	long long x, y, i, p;
	real weight;
	std::map<std::string, long long> dict;
	n_vertices = 0;
	FILE *fin = fopen(infile, "rb");
	if (fin == NULL)
	{
		printf("\nFile not found!\n");
		return;
	}
	printf("Reading input file %s ......%c", infile, 13);
	while (fscanf(fin, "%s%s%f", w1, w2, &weight) == 3)
	{
		if (!dict.count(w1)) { dict[w1] = n_vertices++; names.push_back(w1); }
		if (!dict.count(w2)) { dict[w2] = n_vertices++; names.push_back(w2); }
		x = dict[w1];
		y = dict[w2];
		edge_from.push_back(x);
		edge_to.push_back(y);
		edge_weight.push_back(weight);
		next.push_back(-1);
		++n_edge;
		if (n_edge % 5000 == 0)
		{
			printf("Reading input file %s ...... %lldK edges%c", infile, n_edge / 1000, 13);
			fflush(stdout);
		}
	}
	fclose(fin);
	delete[] w1;
	delete[] w2;

	head = new long long[n_vertices];
	for (i = 0; i < n_vertices; ++i) head[i] = -1;
	for (p = 0; p < n_edge; ++p)
	{
		x = edge_from[p];
		next[p] = head[x];
		head[x] = p;
	}
	printf("\nTotal vertices : %lld\tTotal edges : %lld\n", n_vertices, n_edge);
}

void LargeVis::save(char *outfile,char *outfile1) {
	FILE *fout = fopen(outfile, "wb");

	long long NxD = n_vertices * out_dim;
	fprintf(fout, "%lld %lld\n", n_vertices, out_dim);
	for (long long map = 0; map < no_maps; ++map) {
		for (long long i = 0; i < n_vertices; ++i) {
			if (names.size()) fprintf(fout, "%s ", names[i].c_str());
			for (long long j = 0; j < out_dim; ++j) {
				if (j) fprintf(fout, " ");
				fprintf(fout, "%.6f", vis[i * out_dim + j + NxD * map]);
			}
			fprintf(fout, "\n");
		}
	}
	fclose(fout);
	FILE *fout1 = fopen(outfile1,"wb");
	for (long long i = 0; i < n_vertices; ++i) {
		for (int map = 0; map < no_maps; ++map) {
			if(map) fprintf(fout1," ");
			fprintf(fout1,"%.6f",weights[i * no_maps + map]);
		}

		fprintf(fout1,"\n");
	}
	fclose(fout1);
}

real *LargeVis::get_ans()
{
	return vis;
}

long long LargeVis::get_n_vertices()
{
	return n_vertices;
}

long long LargeVis::get_out_dim()
{
	return out_dim;
}

void LargeVis::normalize()
{
    	printf("Normalizing ......"); fflush(stdout);
	real *mean = new real[n_dim];
	for (long long i = 0; i < n_dim; ++i) mean[i] = 0;
	for (long long i = 0, ll = 0; i < n_vertices; ++i, ll += n_dim)
	{
		for (long long j = 0; j < n_dim; ++j)
			mean[j] += vec[ll + j];
	}
	for (long long j = 0; j < n_dim; ++j)
		mean[j] /= n_vertices;
	real mX = 0;
	for (long long i = 0, ll = 0; i < n_vertices; ++i, ll += n_dim)
	{
		for (long long j = 0; j < n_dim; ++j)
		{
			vec[ll + j] -= mean[j];
			if (fabs(vec[ll + j]) > mX)	mX = fabs(vec[ll + j]);
		}
	}
	for (long long i = 0; i < n_vertices * n_dim; ++i)
		vec[i] /= mX;
	delete[] mean;
	printf(" Done.\n");
}

real LargeVis::CalcDist(long long x, long long y)
{
	real ret = 0;
	long long i, lx = x * n_dim, ly = y * n_dim;
	for (i = 0; i < n_dim; ++i)
		ret += (vec[lx + i] - vec[ly + i]) * (vec[lx + i] - vec[ly + i]);
	return ret;
}
real LargeVis::CalcDistVis(long long x, long long y) {
	real ret = 0;
	long long i, lx = x * out_dim, ly = y * out_dim;
	for (i = 0; i < out_dim; ++i)
		ret += (vis[lx + i] - vis[ly + i]) * (vis[lx + i] - vis[ly + i]);
	return ret;
}

void LargeVis::init_weight_table() {
	real value = 1.0 / no_maps;
	printf("n_smapls is %lld\n", n_samples);
	for (long long map = 0;map < no_maps; ++map) {
		for (long long i = 0; i < n_vertices; ++i) {
			weights[i * no_maps + map] = value;
		}
	}
}
//根据 P 矩阵的值构建边采样的概率
void LargeVis::init_alias_table()
{
	alias = new long long[n_edge];
	prob = new real[n_edge];

	real *norm_prob = new real[n_edge];
	long long *large_block = new long long[n_edge];
	long long *small_block = new long long[n_edge];

	real sum = 0;
	long long cur_small_block, cur_large_block;
	long long num_small_block = 0, num_large_block = 0;

	for (long long k = 0; k < n_edge; ++k) sum += edge_weight[k];
	for (long long k = 0; k < n_edge; ++k) norm_prob[k] = edge_weight[k] * n_edge / sum;

	for (long long k = n_edge - 1; k >= 0; --k)
	{
		if (norm_prob[k] < 1)
			small_block[num_small_block++] = k;
		else
			large_block[num_large_block++] = k;
	}

	while (num_small_block && num_large_block)
	{
		cur_small_block = small_block[--num_small_block];
		cur_large_block = large_block[--num_large_block];
		prob[cur_small_block] = norm_prob[cur_small_block];
		alias[cur_small_block] = cur_large_block;
		norm_prob[cur_large_block] = norm_prob[cur_large_block] + norm_prob[cur_small_block] - 1;
		if (norm_prob[cur_large_block] < 1)
			small_block[num_small_block++] = cur_large_block;
		else
			large_block[num_large_block++] = cur_large_block;
	}

	while (num_large_block) prob[large_block[--num_large_block]] = 1;
	while (num_small_block) prob[small_block[--num_small_block]] = 1;

	delete[] norm_prob;
	delete[] small_block;
	delete[] large_block;
}

long long LargeVis::sample_an_edge(real rand_value1, real rand_value2)
{
	long long k = (long long)((n_edge - 0.1) * rand_value1);
	return rand_value2 <= prob[k] ? k : alias[k];
}

void LargeVis::annoy_thread(int id)
{
	long long lo = id * n_vertices / n_threads;
	long long hi = (id + 1) * n_vertices / n_threads;
	AnnoyIndex<int, real, Euclidean, Kiss64Random> *cur_annoy_index = NULL;
	if (id > 0)
	{
		cur_annoy_index = new AnnoyIndex<int, real, Euclidean, Kiss64Random>(n_dim);
		cur_annoy_index->load("annoy_index_file");
	}
	else
		cur_annoy_index = annoy_index;
	for (long long i = lo; i < hi; ++i)
	{
		cur_annoy_index->get_nns_by_item(i, n_neighbors + 1, (n_neighbors + 1) * n_trees, &knn_vec[i], NULL);
		for (long long j = 0; j < knn_vec[i].size(); ++j)
			if (knn_vec[i][j] == i)
			{
				knn_vec[i].erase(knn_vec[i].begin() + j);
				break;
			}
	}
	if (id > 0) delete cur_annoy_index;
}

void *LargeVis::annoy_thread_caller(void *arg)
{
	LargeVis *ptr = (LargeVis*)(((arg_struct*)arg)->ptr);
	ptr->annoy_thread(((arg_struct*)arg)->id);
	pthread_exit(NULL);
}

void LargeVis::run_annoy()
{
    	printf("Running ANNOY ......"); fflush(stdout);
	annoy_index = new AnnoyIndex<int, real, Euclidean, Kiss64Random>(n_dim);
	for (long long i = 0; i < n_vertices; ++i)
		annoy_index->add_item(i, &vec[i * n_dim]);
	annoy_index->build(n_trees);
	if (n_threads > 1) annoy_index->save("annoy_index_file");
	knn_vec = new std::vector<int>[n_vertices];

	pthread_t *pt = new pthread_t[n_threads];
	for (int j = 0; j < n_threads; ++j) pthread_create(&pt[j], NULL, LargeVis::annoy_thread_caller, new arg_struct(this, j));
	for (int j = 0; j < n_threads; ++j) pthread_join(pt[j], NULL);
	delete[] pt;
    	delete annoy_index; annoy_index = NULL;
	printf(" Done.\n");
}

void LargeVis::propagation_thread(int id)
{
	long long lo = id * n_vertices / n_threads;
	long long hi = (id + 1) * n_vertices / n_threads;
	int *check = new int[n_vertices];
	std::priority_queue< pair<real, int> > heap;
	long long x, y, i, j, l1, l2;
	for (x = 0; x < n_vertices; ++x) check[x] = -1;
	for (x = lo; x < hi; ++x)
	{
		check[x] = x;
		std::vector<int> &v1 = old_knn_vec[x];
		l1 = v1.size();
		for (i = 0; i < l1; ++i)
		{
			y = v1[i];
			check[y] = x;
			heap.push(std::make_pair(CalcDist(x, y), y));
			if (heap.size() == n_neighbors + 1) heap.pop();
		}
		for (i = 0; i < l1; ++i)
		{
			std::vector<int> &v2 = old_knn_vec[v1[i]];
			l2 = v2.size();
			for (j = 0; j < l2; ++j) if (check[y = v2[j]] != x)
			{
				check[y] = x;
				heap.push(std::make_pair(CalcDist(x, y), (int)y));
				if (heap.size() == n_neighbors + 1) heap.pop();
			}
		}
		while (!heap.empty())
		{
			knn_vec[x].push_back(heap.top().second);
			heap.pop();
		}
	}
	delete[] check;
}

void *LargeVis::propagation_thread_caller(void *arg)
{
	LargeVis *ptr = (LargeVis*)(((arg_struct*)arg)->ptr);
	ptr->propagation_thread(((arg_struct*)arg)->id);
	pthread_exit(NULL);
}

void LargeVis::run_propagation()
{
	for (int i = 0; i < n_propagations; ++i)
	{
		printf("Running propagation %d/%d%c", i + 1, n_propagations, 13);
		fflush(stdout);
		old_knn_vec = knn_vec;
		knn_vec = new std::vector<int>[n_vertices];
		pthread_t *pt = new pthread_t[n_threads];
		for (int j = 0; j < n_threads; ++j) pthread_create(&pt[j], NULL, LargeVis::propagation_thread_caller, new arg_struct(this, j));
		for (int j = 0; j < n_threads; ++j) pthread_join(pt[j], NULL);
		delete[] pt;
		delete[] old_knn_vec;
		old_knn_vec = NULL;
	}
	printf("\n");
}

// 计算 P 概率矩阵 => 存储在 edgeweight中 => 一个线程只存一点?
void LargeVis::compute_similarity_thread(int id)
{
	long long lo = id * n_vertices / n_threads;
	long long hi = (id + 1) * n_vertices / n_threads;
	long long x, iter, p;
	real beta, lo_beta, hi_beta, sum_weight, H, tmp;
	for (x = lo; x < hi; ++x)
	{
		beta = 1;
		lo_beta = hi_beta = -1;
		for (iter = 0; iter < 200; ++iter)
		{
			H = 0;
            		sum_weight = FLT_MIN;
			for (p = head[x]; p >= 0; p = next[p])
			{
				sum_weight += tmp = exp(-beta * edge_weight[p]);
				H += beta * (edge_weight[p] * tmp);
			}
			H = (H / sum_weight) + log(sum_weight);
			if (fabs(H - log(perplexity)) < 1e-5) break;
			if (H > log(perplexity))
			{
				lo_beta = beta;
				if (hi_beta < 0) beta *= 2; else beta = (beta + hi_beta) / 2;
			}
			else{
				hi_beta = beta;
				if (lo_beta < 0) beta /= 2; else beta = (lo_beta + beta) / 2;
			}
            		if(beta > FLT_MAX) beta = FLT_MAX;
        	}
		for (p = head[x], sum_weight = FLT_MIN; p >= 0; p = next[p])
		{
			sum_weight += edge_weight[p] = exp(-beta * edge_weight[p]);
		}
		for (p = head[x]; p >= 0; p = next[p])
		{
			edge_weight[p] /= sum_weight;
		}
	}
}

void *LargeVis::compute_similarity_thread_caller(void *arg)
{
	LargeVis *ptr = (LargeVis*)(((arg_struct*)arg)->ptr);
	ptr->compute_similarity_thread(((arg_struct*)arg)->id);
	pthread_exit(NULL);
}

void LargeVis::search_reverse_thread(int id)
{
	long long lo = id * n_vertices / n_threads;
	long long hi = (id + 1) * n_vertices / n_threads;
	long long x, y, p, q;
	for (x = lo; x < hi; ++x)
	{
		for (p = head[x]; p >= 0; p = next[p])
		{
			y = edge_to[p];
			for (q = head[y]; q >= 0; q = next[q])
			{
				if (edge_to[q] == x) break;
			}
			reverse[p] = q;
		}
	}
}

void *LargeVis::search_reverse_thread_caller(void *arg)
{
	LargeVis *ptr = (LargeVis*)(((arg_struct*)arg)->ptr);
	ptr->search_reverse_thread(((arg_struct*)arg)->id);
	pthread_exit(NULL);
}

void LargeVis::compute_similarity()
{
    	printf("Computing similarities ......"); fflush(stdout);
	n_edge = 0;
	head = new long long[n_vertices];
	long long i, x, y, p, q;
	real sum_weight = 0;
	for (i = 0; i < n_vertices; ++i) head[i] = -1;
	for (x = 0; x < n_vertices; ++x)
	{
		for (i = 0; i < knn_vec[x].size(); ++i)
		{
			edge_from.push_back((int)x);
			edge_to.push_back((int)(y = knn_vec[x][i]));
			edge_weight.push_back(CalcDist(x, y));
			next.push_back(head[x]);
			reverse.push_back(-1);
			head[x] = n_edge++;
		}
	}
		//for visualize result
    	delete[] vec; vec = NULL;
    	delete[] knn_vec; knn_vec = NULL;
	pthread_t *pt = new pthread_t[n_threads];
	for (int j = 0; j < n_threads; ++j) pthread_create(&pt[j], NULL, LargeVis::compute_similarity_thread_caller, new arg_struct(this, j));
	for (int j = 0; j < n_threads; ++j) pthread_join(pt[j], NULL);
	delete[] pt;

	pt = new pthread_t[n_threads];
	for (int j = 0; j < n_threads; ++j) pthread_create(&pt[j], NULL, LargeVis::search_reverse_thread_caller, new arg_struct(this, j));
	for (int j = 0; j < n_threads; ++j) pthread_join(pt[j], NULL);
	delete[] pt;

	for (x = 0; x < n_vertices; ++x)
	{
		for (p = head[x]; p >= 0; p = next[p])
		{
			y = edge_to[p];
			q = reverse[p];
			if (q == -1)
			{
				edge_from.push_back((int)y);
				edge_to.push_back((int)x);
				edge_weight.push_back(0);
				next.push_back(head[y]);
				reverse.push_back(p);
				q = reverse[p] = head[y] = n_edge++;
			}
			if (x > y){
				sum_weight += edge_weight[p] + edge_weight[q];
				edge_weight[p] = edge_weight[q] = (edge_weight[p] + edge_weight[q]) / 2;
			}
		}
	}
	printf(" Done.\n");
}

void LargeVis::test_accuracy()
{
	long long test_case = 100;
	std::priority_queue< pair<real, int> > *heap = new std::priority_queue< pair<real, int> >;
	long long hit_case = 0, i, j, x, y;
	for (i = 0; i < test_case; ++i)
	{
		x = floor(gsl_rng_uniform(gsl_r) * (n_vertices - 0.1));
		for (y = 0; y < n_vertices; ++y) if (x != y)
		{
			heap->push(std::make_pair(CalcDist(x, y), y));
			if (heap->size() == n_neighbors + 1) heap->pop();
		}
		while (!heap->empty())
		{
			y = heap->top().second;
			heap->pop();
			for (j = 0; j < knn_vec[x].size(); ++j) if (knn_vec[x][j] == y)
				++hit_case;
		}
	}
    	delete heap;
	printf("Test knn accuracy : %.2f%%\n", hit_case * 100.0 / (test_case * n_neighbors));
}

void LargeVis::computeSquaredEuclideanDistance(real *vec1, long long N, long long D, real *DD) {
	const real * XnD = vec1;
	for(int n = 0; n < N; ++n, XnD += D) {
		const real * XmD = XnD + D;
		real * curr_elem = &DD[n*N + n];
		*curr_elem = 0.0;
		real * curr_elem_sym = curr_elem + N;
		for(int m = n + 1; m < N; ++m, XmD+=D, curr_elem_sym+=N) {
			*(++curr_elem) = 0.0;
			for(int d = 0; d < D; ++d) {
				*curr_elem += (XnD[d] - XmD[d]) * (XnD[d] - XmD[d]);
			}
			*curr_elem_sym = *curr_elem;
		}
	}
}
void LargeVis::computeSquaredEuclideanDistanceNo_maps(real * X, int N, int D, real * DD, int no_maps_) {
	if (no_maps_ == 1) {
		computeSquaredEuclideanDistance(X,N, D, DD);
		return;
	}
	for (int map = 0; map < no_maps_; ++map) {
		const real* XnD = X + map * N * D;
		for (int n = 0; n < N; ++n, XnD += D) {
			const real *XmD = XnD + D;
			real *curr_elem = &DD[n * N + n + map * N * N];
			*curr_elem = 0.0;
			real *curr_elem_sym = curr_elem + N;
			for (int m = n + 1; m < N; ++m, XmD += D, curr_elem_sym += N) {
				*(++curr_elem) = 0.0;
				for (int d = 0; d < D; ++d) {
					*curr_elem += (XnD[d] - XmD[d]) * (XnD[d] - XmD[d]);
				}
				*curr_elem_sym = *curr_elem;
			}
		}
	}
}
void LargeVis::computeGaussianPerplexity(real* P, real perplexity) {
	// Compute the squared Euclidean distance matrix
	real * DD = (real *) malloc(n_vertices * n_vertices * sizeof(real));
	if(DD == NULL) { printf("Memory allocation failed!\n"); exit(1); }
	load_from_file("/Users/bigheiniu/course/graduate_pro/tsne/spark-tsne-master/data/mnist/mnistless.txt");
	computeSquaredEuclideanDistance(vec, n_vertices, n_dim, DD);

	// Compute the Gaussian kernel row by row
	int nN = 0;
	for(int n = 0; n < n_vertices; n++) {

		// Initialize some variables
		bool found = false;
		real beta = 1.0;
		real min_beta = -DBL_MAX;
		real max_beta =  DBL_MAX;
		real tol = 1e-5;
		real sum_P;

		// Iterate until we found a good perplexity
		int iter = 0;
		while(!found && iter < 200) {

			// Compute Gaussian kernel row
			for(int m = 0; m < n_vertices; m++) P[nN + m] = exp(-beta * DD[nN + m]);
			P[nN + n] = DBL_MIN;

			// Compute entropy of current row
			sum_P = DBL_MIN;
			for(int m = 0; m < n_vertices; m++) sum_P += P[nN + m];
			real H = 0.0;
			for(int m = 0; m < n_vertices; m++) H += beta * (DD[nN + m] * P[nN + m]);
			H = (H / sum_P) + log(sum_P);

			// Evaluate whether the entropy is within the tolerance level
			real Hdiff = H - log(perplexity);
			if(Hdiff < tol && -Hdiff < tol) {
				found = true;
			}
			else {
				if(Hdiff > 0) {
					min_beta = beta;
					if(max_beta == DBL_MAX || max_beta == -DBL_MAX)
						beta *= 2.0;
					else
						beta = (beta + max_beta) / 2.0;
				}
				else {
					max_beta = beta;
					if(min_beta == -DBL_MAX || min_beta == DBL_MAX)
						beta /= 2.0;
					else
						beta = (beta + min_beta) / 2.0;
				}
			}

			// Update iteration counter
			iter++;
		}

		// Row normalize P
		for(int m = 0; m < n_vertices; m++) P[nN + m] /= sum_P;
		nN += n_vertices;
	}

	// Clean up memory
	free(DD); DD = NULL;

}
void LargeVis::computeQ(real *Q) {
	if (Q == NULL) {
		printf("not enough memory for Q\n");
		exit(-1);
	} else {
		for (long long i = 0; i < n_vertices * n_vertices; ++i) {
			Q[i] = 0;
		}
	}
	real *DD = (real*)malloc(n_vertices * n_vertices * no_maps * sizeof(real));
	real * proportions = (real *) malloc(n_vertices * no_maps * sizeof(real));
	if (DD == NULL || proportions == NULL) {
		printf("not enough memory in computeQ \n");
		exit(-1);
	}
	computeSquaredEuclideanDistanceNo_maps(vis,n_vertices,out_dim,DD,no_maps);

	real sum_proportion,sum_Q=0;
	long long weightIndex = 0,nN = 0;
	for (long long n = 0; n < n_vertices; ++n) {
		sum_proportion = 0;

		for (int map = 0; map < no_maps; ++map) {
			sum_proportion += expf(weights[weightIndex + map]);
		}
		for (int map = 0; map < no_maps; ++map) {
			proportions[weightIndex + map] = expf(weights[weightIndex + map]) / sum_proportion;
		}
		weightIndex += no_maps;
	}

	for(long long n = 0; n < n_vertices; n++) {
		for(long long m = 0; m < n_vertices; m++) {
			for( int map = 0; map < no_maps; ++map) {
				if (n != m) Q[nN + m] += proportions[n*no_maps+map] * proportions[m*no_maps+map] / (1 + DD[nN + m + map * n_vertices * n_vertices]);
			}
			if (n!=m) {
				sum_Q += Q[nN + m];
			}
			else Q[nN + m] = DBL_MIN;
		}
		nN += n_vertices;
	}
	for(int i = 0; i < n_vertices * n_vertices; i++) Q[i] /= sum_Q;


	free(DD);DD=NULL;
	free(proportions);proportions=NULL;
}
void LargeVis::computeP(real * P) {
    if (P == NULL) {
        printf("not enough memory for P\n");
        exit(-1);
    } else {
        for (long long i = 0; i < n_vertices * n_vertices; ++i) {
            P[i] = 0;
        }
    }
	computeGaussianPerplexity(P, perplexity);
    int nN = 0;
    for(int n = 0; n < n_vertices; n++) {
        int mN = (n + 1) * n_vertices;
        for(int m = n + 1; m < n_vertices; m++) {
            P[nN + m] += P[mN + n];
            P[mN + n]  = P[nN + m];
            mN += n_vertices;
        }
        nN += n_vertices;
    }
    real sum_P = 0;
    for(int i = 0; i < n_vertices * n_vertices; i++) sum_P += P[i];
    for(int i = 0; i < n_vertices * n_vertices; i++) P[i] /= sum_P;
}

void LargeVis::test_visualize() {
	real *P = (real *)malloc(n_vertices * n_vertices * sizeof(real));
	computeP(P);
	real *Q = (real *) malloc(n_vertices * n_vertices * sizeof(real));
	computeQ(Q);
	int count = 0;
	double value1 ;
	int index1;
	double value2 ;
	int index2;
	for(int i = 0; i < n_vertices; ++i) {
		value1 = -1;
		value2 = -1;
		index1 = -1;
		index2 = -1;
		for (int j = 0; j < n_vertices; ++j) {
			if (value1 < P[i * n_vertices + j]) {
				value1 = P[i * n_vertices + j];
				index1 = j;
			}
			if (value2 < Q[i * n_vertices + j]) {
				value2 = Q[i * n_vertices + j];
				index2 = j;
			}
		}
		if (index1 == index2) {
			count++;
		}
	}
	printf("NPR value is %f\n",(count * 1.0 / n_vertices));
	free(P); P = NULL;
	free(Q); Q = NULL;
}

void LargeVis::construt_knn()
{
	normalize();
	run_annoy();
	run_propagation();
	test_accuracy();
	compute_similarity();

	/*FILE *fout = fopen("knn_graph.txt", "wb");
	for (long long p = 0; p < n_edge; ++p)
	{
		fprintf(fout, "%lld %lld ", edge_from[p], edge_to[p]);
		real tmp = edge_weight[p];
		fwrite(&tmp, sizeof(real), 1, fout);
		fprintf(fout, "\n");
	}
	fclose(fout);*/
}

//负采样边的构建
void LargeVis::init_neg_table()
{
	long long x, p, i;
	neg_size = 1e8;
	reverse.clear(); vector<long long> (reverse).swap(reverse);
	// dd 应该是度
	real sum_weights = 0, dd, *real1 = new real[n_vertices];
	for (i = 0; i < n_vertices; ++i) real1[i] = 0;
	for (x = 0; x < n_vertices; ++x)
	{
		for (p = head[x]; p >= 0; p = next[p])
		{
			real1[x] += edge_weight[p];
		}
		sum_weights += real1[x] = pow(real1[x], 0.75);
	}
    	next.clear(); vector<long long> (next).swap(next);
    	delete[] head; head = NULL;
	neg_table = new int[neg_size];
	dd = real1[0];
	for (i = x = 0; i < neg_size; ++i)
	{
		neg_table[i] = x;
		if (i / (real)neg_size > dd / sum_weights && x < n_vertices - 1)
		{
			dd += real1[++x];
		}
	}
	delete[] real1;
}

//check ok!
real LargeVis::proportions_sum(long long index) {
    real sum = 0;
    long long mapIndex = index * no_maps;
    for (long long map = 0; map < no_maps; ++map) {
    	sum += expf(-weights[mapIndex + map]);
    }
    if (isnan(sum)||isinf(sum)) {
		for (long long map = 0; map < no_maps; ++map) {
			printf("weight %lld map NO.%lld is %lf\n",index,map,-weights[index * no_maps + map]);
            fflush(stdout);
		}
    }
    return sum;
}

// dOdW = dOdPi * dPidW, 在这里提前计算dPdW
void LargeVis::Calc_pi_weight(real *piArray,real *dPdWArray, long long index) {
    real sumValue = 0.0001;
    long long flagIndex = index * no_maps;
    for(long long map = 0; map < no_maps; ++map) {
        piArray[map] = exp(weights[flagIndex + map]);
		sumValue += piArray[map];
    }
    for(long long map = 0; map < no_maps; ++map) {
	 dPdWArray[map] = (sumValue - piArray[map]) * piArray[map] / (sumValue * sumValue);
	if(dPdWArray[map] == 0) {
		for(int i = 0; i < no_maps; ++i) {
			printf("dpdw is zero, the map NO.%d is %f and index is %lld\n",i,weights[flagIndex+ i],(flagIndex/no_maps));	
	}
	}
    	piArray[map] /= sumValue;
    }
}

//compute t and a, b;
//check ok;
real LargeVis::Calc_T_A_B(long long indexX, long long indexY,real *pi_i , real *pi_j,real* a, real* b) {
    real tij = 0;
    long long NxD = n_vertices * out_dim;
    long long mapIndex = 0;
    int d;
    real pi_ij;
	real distance;
    for(long long map = 0; map < no_maps; ++map, mapIndex += NxD) {
        //weights 太小导致溢出;
     	pi_ij = pi_i[map] * pi_j[map];
        distance = 0;
        for(d = 0; d < out_dim; ++d){
        	// 计算 ( y_i(m) - y_j(m) )^2
			distance += (vis[mapIndex + indexX + d] - vis[mapIndex + indexY + d]) * (vis[mapIndex+indexX + d] - vis[mapIndex+indexY + d]); }
        b[map] = pi_ij / (1.0 + distance);
        a[map]= b[map] / (1.0 + distance);
		tij += b[map];

    }
    if((tij- 0) < 0.000000000001) {

			printf("map NO.%lld tij is %f, value is %f, pi is %f, indexX is %lld, indexY is %lld\n",4,tij,distance,pi_ij,indexX,indexY);
			printf("expindexX value is %f, expindexY value is %f\n",(weights[indexX/out_dim *no_maps]),(weights[indexY/out_dim*no_maps]));
		}
    return tij;
}


void LargeVis::visualize_thread(int id)
{
	//TODO: set dCdW clip
	printf("thread id: %d\n",id);
	fflush(stdout);
	long long edge_count = 0, last_edge_count = 0;
	long long x, y, p, lx, ly, i, j, map,d,mapIndex;
	long long DimxMap = out_dim * no_maps;
	real g, gg, cur_alpha = initial_alpha;
	real dCdW = 0.0;
	real *sumEdgeProportions = new real[no_maps];
	real *A = new real[no_maps];
	real *B = new real[no_maps];
	real *edgeProportions = new real[no_maps];
	real *cur = new real[out_dim];
	real *err = new real[DimxMap];
	real weight_clip = 10;
	//设置这两个变量还是存在不同步的问题;
	real *pi_x = new real[no_maps];
    real *pi_y = new real[no_maps];
    real *dPdW_x = new real[no_maps];
    real *dPdW_y = new real[no_maps];

	if(A ==NULL || B == NULL || edgeProportions == NULL || sumEdgeProportions == NULL || cur == NULL || err == NULL) {
		printf("not enough memory \n");
		exit(-1);
	}
	long long NxD = n_vertices * out_dim;
	// Puzzle:grad_clip 是什么作用
	real grad_clip = 5.0;
	while (1) {
        // 在 mm tsne 中边采样采的是高维矩阵的边, 所以跟图的数量没有关系
        if (edge_count > (n_samples / n_threads + 2)) {
        	printf("edge count is %d\n",edge_count);
        	printf("value is %lld\n",n_samples/n_threads);
        	break;
        }
        if (edge_count - last_edge_count > 10000) {
            // 根据迭代对模型进行放缩
			edge_count_actual += edge_count - last_edge_count;
            last_edge_count = edge_count;
            cur_alpha = initial_alpha * (1 - edge_count_actual / (n_samples + 1.0));
            if (cur_alpha < initial_alpha * 0.0001) cur_alpha = initial_alpha * 0.0001;
            printf("%cFitting model\tAlpha: %f Progress: %.3lf%% edge count%lld", 13, cur_alpha,
                   (real) edge_count_actual / (real) (n_samples + 1) * 100,edge_count_actual);
			fflush(stdout);
        }
        p = sample_an_edge(gsl_rng_uniform(gsl_r), gsl_rng_uniform(gsl_r));
        x = edge_from[p];
	if(x > n_vertices) {
	   printf("sample error\n");
	  exit(-1);
	}
		Calc_pi_weight(pi_x, dPdW_x,x);
		lx = x * out_dim;
		y = edge_to[p];

        for (map = 0; map < no_maps; ++map){
        	//存储dCdP;
        	sumEdgeProportions[map] = 0;
        }

        for(i = 0; i < DimxMap; ++i) {
        	err[i] = 0;
        }

		for (i = 0; i < n_negatives + 1; ++i) {
			if (i > 0) {
				y = neg_table[(unsigned long long) floor(gsl_rng_uniform(gsl_r) * (neg_size - 0.1))];
				if (y == edge_to[p] || y == edge_from[p]) continue;
			}
			ly = y * out_dim;
			if(lx > n_vertices * out_dim) {
			  printf("dimention error %lld and x is %lld\n",lx,x);
			  printf("vertices is %lld\n",n_vertices);
			  printf("out_dim is %lld\n",out_dim);
			 exit(-1);
			}
			Calc_pi_weight(pi_y, dPdW_y, y);
			//sum 为 求导很大那一坨东西
					// A = π_i^m * π_i^m / (1+d_ij^m)^2
			real T =Calc_T_A_B(lx,ly, pi_x, pi_y, A, B);
			if(isnan(T)) {
				printf("T is nan\n");
                fflush(stdout);
				exit(-1);
			}
			// dCdW 公式计算需要对 dCdP(i,m)进行求和

			for (map = 0,mapIndex=0; map < no_maps; ++map,mapIndex += NxD) {
				for (d = 0; d < out_dim; ++d) {
					cur[d] = vis[lx+mapIndex+d];
				}
				if(i == 0) {
					g = -2 / T * A[map];
					if (isnan(g)) {
						printf("need1 soft and T is %f\n",T);
						fflush(stdout);
						exit(-1);
					}
					// dCdP 是半成品, 在 mm-tSNE中, 它为π_i^m' * π_j^m' / (1 + d_ij^m')
					edgeProportions[map] = 1 / T * B[map];
				}
				else {
					g = 2 * gamma / (1 - T) * A[map];
					if (isnan(g)) {
						printf("need soft and T is %f\n",T);
                        fflush(stdout);
						exit(-1);
					}
					// dCdP 是半成品, 在 mm-tSNE中, 它为π_i^m' * dCdP_i^m'
					edgeProportions[map] = -1 * gamma / (1 - T) * B[map] ;
				}
				if(edgeProportions[map] > 1.6) edgeProportions[map] = 1.6;
				if(edgeProportions[map] < -1.6) edgeProportions[map] = -1.6;
				sumEdgeProportions[map] += edgeProportions[map];
				// 计算 dCdW_y^m 需要
				for(j = 0; j < out_dim; ++j)
					gg = g * (cur[j] - vis[ly +mapIndex + j]);
					if (gg > grad_clip) gg = grad_clip;
					if (gg < -grad_clip) gg = -grad_clip;
					err[j + map * out_dim] += gg * cur_alpha;
					if (isnan(err[j + map * out_dim])) {
						printf("visX is nan, g is %f, cur[%lld] is %f, vis[%lld] is %f\n",g, j, cur[j],(ly+j),vis[ly+j]);
					}

					gg = g * (vis[ly +mapIndex+ j] - cur[j]);
					if (gg > grad_clip) gg = grad_clip;
					if (gg < -grad_clip) gg = -grad_clip;
					vis[ly +mapIndex+ j] += gg * cur_alpha;
					weights[y * no_maps + map] += edgeProportions[map] / pi_y[map] * dPdW_y[map];
					if(weights[y * no_maps + map] > weight_clip) {
						weights[y * no_maps + map] = weight_clip;
					}
					if(weights[y * no_maps + map] < -weight_clip) {
						weights[y * no_maps + map] = -weight_clip;
					}
					if(pi_y[map] == 0) {
					 printf("pi_y map NO.%d is 0\n",map);
					exit(-1);	
					}
					if(dPdW_y[map]== 0) {
					printf("dPdW_y NO.%d is empty\n",map);
					exit(-1);
					}
					if (isnan(vis[ly+mapIndex+j])) {
						printf("errorerror %lld; edge is %lld\n",ly+j,edge_count);
						exit(-1);
					}
				}
			}
		}
		for (map = 0,mapIndex = 0; map < no_maps; ++map, mapIndex += NxD) {
			for (j = 0; j < out_dim; ++j) {
			    vis[x * out_dim + mapIndex + j] += err[out_dim * map + j];
                if (isnan(vis[x * out_dim + mapIndex + j])) {
                    printf("err NO.%d dimention NO.%d is %f\n",map,j,err[out_dim*map +j]);
                }
			}
            weights[x * no_maps + map] += sumEdgeProportions[map] / pi_x[map] *dPdW_x[map];
		if(weights[x * no_maps + map]  > weight_clip) weights[x * no_maps + map] = weight_clip;
			if(weights[x * no_maps + map]  < -weight_clip) weights[x * no_maps + map] = -weight_clip;
		if(isnan(weights[x * no_maps + map])) {
			printf("weights nan, x %lld, map %d\n",x,map);
		}
		if( std::abs(sumEdgeProportions[map] / pi_x[map] * dPdW_x[map])< 0.000001) {
			printf("edge:%lld  gradient error\n",edge_count);
			printf("value is %f\n",(sumEdgeProportions[map] / pi_x[map] * dPdW_x[map]));
			printf("pi_x value %f, sumEdgeProportions value %f, dPdW_x value %f\n",pi_x[map],sumEdgeProportions[map],dPdW_x[map]);
		}
		}
        ++edge_count;
    }
	delete[] cur;
	delete[] err;
	delete[] A;
	delete[] B;
	delete[] pi_x;
	delete[] pi_y;
	delete[] dPdW_x;
	delete[] dPdW_y;
	delete[] sumEdgeProportions;
	delete[] edgeProportions;

}
void *LargeVis::visualize_thread_caller(void *arg)
{
	LargeVis *ptr = (LargeVis*)(((arg_struct*)arg)->ptr);
	ptr->visualize_thread(((arg_struct*)arg)->id);
	pthread_exit(NULL);
}

void LargeVis::visualize()

{
//	printf("fuck! %lld\n",n_samples);
	long long i;
	vis = new real[n_vertices * out_dim * no_maps];
	weights = new real[n_vertices * no_maps];
	for (i = 0; i < n_vertices * out_dim * no_maps; ++i) vis[i] = (gsl_rng_uniform(gsl_r) - 0.5) / out_dim * 0.0001;
	for (i = 0;  i< no_maps * n_vertices; ++i) weights[i] = ((gsl_rng_uniform(gsl_r) - 0.5 ) / no_maps) * 0.0001;
	init_neg_table();
	init_alias_table();
	edge_count_actual = 0;
	//test set thread to 1;
		n_threads = 1;
	pthread_t *pt = new pthread_t[n_threads];
	for (int j = 0; j < n_threads; ++j) pthread_create(&pt[j], NULL, LargeVis::visualize_thread_caller, new arg_struct(this, j));
	for (int j = 0; j < n_threads; ++j) pthread_join(pt[j], NULL);
	delete[] pt;
//	printf("\n");
}

void LargeVis::run(long long out_d, long long n_thre, long long n_samp, long long n_prop, real alph, long long n_tree, long long n_nega, long long n_neig, real gamm, real perp,long long no_m)
{
	clock_t startime,endtime;
	startime = clock();
	gsl_rng_env_setup();
	gsl_T = gsl_rng_rand48;
	gsl_r = gsl_rng_alloc(gsl_T);
	gsl_rng_set(gsl_r, 314159265);

	clean_model();
	if (!vec && !head)
	{
		printf("Missing training data!\n");
		return;
	}
	// maps 为 5
	no_maps = no_m <0 ? 5 : no_m;
//	no_maps = 1;
	out_dim = out_d < 0 ? 2 : out_d;
	initial_alpha = alph < 0 ? 1.0 : alph;
	n_threads = n_thre < 0 ? 8 : n_thre;
	n_samples = n_samp;
	n_negatives = n_nega < 0 ? 5 : n_nega;
	n_neighbors = n_neig < 0 ? 150 : n_neig;
	n_trees = n_tree;
	n_propagations = n_prop < 0 ? 3 : n_prop;
	gamma = gamm < 0 ? 7.0 : gamm;
	perplexity = perp < 0 ? 50.0 : perp;
	if (n_samples < 0)
	{
		if (n_vertices < 10000)
			n_samples = 1000;
		else if (n_vertices < 1000000)
			n_samples = (n_vertices - 10000) * 9000 / (1000000 - 10000) + 1000;
		else n_samples = n_vertices / 100;
	}
	n_samples *= 1000000;
//	n_samples *= 1000;
	if (n_trees < 0)
	{
		if (n_vertices < 100000)
			n_trees = 10;
		else if (n_vertices < 1000000)
			n_trees = 20;
		else if (n_vertices < 5000000)
			n_trees = 50;
		else n_trees = 100;
	}
//	if (vec!=NULL) { clean_graph(); construt_knn(); }
	if (vec) { clean_graph(); construt_knn(); }
	visualize();
	endtime = clock();
	real totaltime =(real)(endtime-startime)/CLOCKS_PER_SEC;
	printf("time for program is %f",totaltime);
}

void LargeVis::load_weight(char *infile) {
   FILE *fin = fopen(infile, "rb");
	if (fin == NULL)
	{
		printf("\nFile not found!\n");
		return;
	}
    	printf("Reading input file %s ......", infile); fflush(stdout);
	fscanf(fin, "%lld%lld", &n_vertices, &no_maps);
	vec = new real[n_vertices * no_maps];
	for (long long i = 0; i < n_vertices; ++i)
	{
		for (long long j = 0; j < no_maps; ++j)
		{
			fscanf(fin, "%f", &weights[i * no_maps + j]);
		}
	}
	fclose(fin);
	printf(" Done.\n");
	printf("Load weight from %s : %lld\tDimension : %lld\n", infile,n_vertices, no_maps);
}
void LargeVis::load_Y(char *infile) {
	FILE *fin = fopen(infile, "rb");
	if (fin == NULL)
	{
		printf("\nFile not found!\n");
		return;
	}
	printf("Reading input file %s ......", infile); fflush(stdout);
	fscanf(fin, "%lld%lld%lld", &n_vertices, &out_dim,&no_maps);
	vis = new real[n_vertices * no_maps * out_dim];
	long long NxD = n_vertices * out_dim;
	for(int map =0; map < no_maps; ++map) {
		for (long long i = 0; i < n_vertices; ++i) {
			for (long long j = 0; j < out_dim; ++j) {
				fscanf(fin, "%f", &vis[i * out_dim + j + map * NxD]);
			}
		}
	}
	fclose(fin);
	printf(" Done.\n");
	printf("Load Y from %s : %lld\tDimension : %lld\n", infile,n_vertices, no_maps);
}

