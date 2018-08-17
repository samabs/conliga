/**
 * Author: Sam Abujudeh
 * 
 * The Assignment Problem Class
 * 
 * Implemented as part of conliga for the purposes of
 * seeking the permutation of labels which minimuses the
 * cost function in the Stephens algorithm.
 * 
 * The Hungarian Algorithm (Munkres) is implemented here
 * and was based on description of the algorithm here
 * http://csclab.murraystate.edu/~bob.pilgrim/445/munkres.html
 * 
 * The Stephens algorithm is implemented in the LabelSwitch Class
 **/


#include "AP.h"

using namespace std;

AP::AP() {}

AP::AP(arma::mat c) {
	cost = c;
	done = false;
	step = 0;
	N = cost.n_rows;
	M = arma::zeros<arma::umat>(N, N);
	rowCover = arma::zeros<arma::ivec>(N);
	colCover = arma::zeros<arma::ivec>(N);
	cpath_0 = 0;
	rpath_0 = 0;
	path.set_size(2 * N, 2);
}

void AP::set(arma::mat c) {
	cost = c;
	done = false;
	step = 0;
	N = cost.n_rows;
	M = arma::zeros<arma::umat>(N, N);
	rowCover = arma::zeros<arma::ivec>(N);
	colCover = arma::zeros<arma::ivec>(N);
	cpath_0 = 0;
	rpath_0 = 0;
	path.set_size(2 * N, 2);
}

// the Munkres assignment algorithm
// Implemented based on description described in:
// http://csclab.murraystate.edu/~bob.pilgrim/445/munkres.html
arma::umat AP::runMunkres() {
	done = false;
	step = 1;
	while (!done) {
		switch (step) {
		
			case 1:
				step_one();
				break;
				
			case 2:
				step_two();
				break;
				
			case 3:
				step_three();
				break;
				
			case 4:
				step_four();
				break;
				
			case 5:
				step_five();
				break;
				
			case 6:
				step_six();
				break;
				
			case 7:
				done = true;
				break;
				
		}
	}
	return M;
}


/**
 * The functions below are used by the Munkres algorithm
 */

void AP::step_one() {
	for (unsigned int r = 0; r < N; ++r) {
		cost.row(r) -= arma::min(cost.row(r));
	}
	step = 2;
}

void AP::step_two() {
	for (unsigned int r = 0; r < N; ++r) {
		for (unsigned int c = 0; c < N; ++c) {
			if (cost.at(r, c) == 0.0 && rowCover.at(r) == 0 && colCover.at(c) == 0) {
				M(r,c) = 1;
				rowCover[r] = 1;
				colCover[c] = 1;
				break; // Only first zero in a row and col
			}
		}
	}

	rowCover.fill(0);
	colCover.fill(0);
	step = 3;
}

void AP::step_three() {
	unsigned int colCount = 0;
	for (unsigned int r = 0; r < N; ++r) {
		for (unsigned int c = 0; c < N; ++c) {
			if (M.at(r,c) == 1) {
				colCover.at(c) = 1;
			}
		}
	}
	for (unsigned int c = 0; c < N; ++c) {
		if (colCover.at(c) == 1) {
			colCount++;
		}
	}
	if (colCount == N) {
		step = 7;
	} else {
		step = 4;
	}
}

void AP::find_noncovered_zero(int &row, int &col) {
	/* used in step four */
	unsigned int r = 0;
	unsigned int c;
	bool complete = false;
	row = -1;
	col = -1;
	while (!complete) {
		c = 0;
		while (true) {
			if (cost.at(r,c) == 0.0 && rowCover.at(r) == 0 && colCover.at(c) == 0) {
				row = r;
				col = c;
				complete = true;
			}
			++c;
			if (c == N || complete) {
				break;
			}
		}
		++r;
		if (r == N) {
			complete = true;
		}
	}
}

bool AP::star_in_row(int &row) {
	/* used in step four */
	bool temp = false;
	for (unsigned int c = 0; c < N; ++c) {
		if (M.at(row,c) == 1) {
			temp = true;
			break;
		}
	}
	return temp;
}

void AP::find_star_in_row(const int &row, int &col) {
	/* used in step four */
	col = -1;
	for (unsigned int c = 0; c < N; ++c) {
		if (M.at(row,c) == 1) {
			col = c;
		}
	}
}

void AP::step_four() {
	int row = -1;
	int col = -1;
	bool complete = false;
	while (!complete) {
		find_noncovered_zero(row, col);
		if (row == -1) {
			complete = true;
			step = 6;
		} else {
			/* uncovered zero */
			M(row, col) = 2;
			if (star_in_row(row)) {
				find_star_in_row(row, col);
				/* Cover the row with the starred zero
                 * and uncover the column with the starred
                 * zero. 
                 */
                 rowCover.at(row) = 1;
                 colCover.at(col) = 0;
			} else {
				/* No starred zero in row with 
                 * uncovered zero 
                 */
                complete = true;
                step = 5;
                rpath_0 = row;
                cpath_0 = col;
			}
		}
	}
}

void AP::find_star_in_col(const int &col, int &row) {
	/* used in step five */ 
    row = -1;
    for (unsigned int r = 0; r < N; ++r) {
        if (M.at(r, col) == 1) {
            row = r;
        }
    }
}

void AP::find_prime_in_row(const int &row, int &col) {
	/* used in step five */ 
    for (unsigned int c = 0; c < N; ++c) {
        if (M.at(row,c) == 2) {
            col = c;
        }
    }
}

void AP::augment_path(const int &path_count) {
	/* used in step five */ 
    for (unsigned int p = 0; p < path_count; ++p) {
        if (M.at(path(p, 0), path(p, 1)) == 1) {
            M.at(path(p, 0), path(p, 1)) = 0;
        } else {
            M.at(path(p, 0), path(p, 1)) = 1;
        }
    }
}

void AP::clear_covers() {
    rowCover.fill(0);
    colCover.fill(0);
}

void AP::erase_primes() {
    for (unsigned int r = 0; r < N; ++r) {
        for (unsigned int c = 0; c < N; ++c) {
            if (M.at(r, c) == 2) {
                M.at(r, c) = 0;
            }
        }
    }
}

void AP::step_five() {
	bool complete = false;
	int row = -1;
	int col = -1;
	unsigned int path_count = 1;
	path.at(path_count - 1, 0) = rpath_0;
	path.at(path_count - 1, 1) = cpath_0;
	while (!complete) {
		find_star_in_col(path.at(path_count - 1, 1), row);
		if (row > -1) {
			/* Starred zero in row 'row' */
            ++path_count;
            path.at(path_count - 1, 0) = row;
            path.at(path_count - 1, 1) = path.at(path_count - 2, 1);
		} else {
			complete = true;
		}
		if (!complete) {
			/* If there is a starred zero find a primed 
             * zero in this row; write index to 'col' */
            find_prime_in_row(path.at(path_count - 1, 0), col);
            ++path_count;
            path.at(path_count - 1, 0) = path.at(path_count - 2, 0);
            path.at(path_count - 1, 1) = col;
		}
	}
	augment_path(path_count);
	clear_covers();
	erase_primes();
	step = 3;
}

void AP::find_smallest(double &minval)
{
    for (unsigned int r = 0; r < N; ++r) {
        for (unsigned int c = 0; c < N; ++c) {
            if (rowCover.at(r) == 0 && colCover.at(c) == 0) {                                                                    
                if (minval > cost.at(r, c)) {
                    minval = cost.at(r, c);
                }
            }
        }
    }
}

void AP::step_six() {
	double minval = DBL_MAX;
    find_smallest(minval);
    for (unsigned int r = 0; r < N; ++r) {
        for (unsigned int c = 0; c < N; ++c) {
            if (rowCover.at(r) == 1) {
                cost.at(r, c) += minval;
            }
            if (colCover.at(c) == 0) {
                cost.at(r, c) -= minval;
            }
        }
    }
    step = 4;
}
