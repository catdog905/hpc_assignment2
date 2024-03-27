#include <iostream>
#include <utility>
#include <vector>
#include <cmath>

using namespace std;

typedef struct Point {
    double x;
    double t;
    double value;
    bool isInitialized;
} point_t;

typedef struct Interval {
    int start, stop;
} interval_t;

// grid[0] is t row for x = 0

vector<vector<point_t>> getGrid(
        interval_t t_domain,
        interval_t x_domain,
        int spatial_step_number,
        int time_step_number) {
    vector<vector<point_t>> grid(
            spatial_step_number + 1,
            vector<point_t>(time_step_number + 1));
    double t_step = (t_domain.stop - t_domain.start) / ((double) time_step_number);
    double x_step = (x_domain.stop - x_domain.start) / ((double) spatial_step_number);
    for (int i = 0; i < spatial_step_number + 1; i++) {
        for (int j = 0; j < time_step_number + 1; j++) {
            grid[i][j] = {i * x_step, j * t_step, 0, false};
        }
    }
    return grid;
};

vector<vector<point_t>> fillGridRowXFixed(vector<vector<point_t>> grid, int row_index, double (*f)(double, double)) {
    vector<vector<point_t>> copy_grid(grid);
    for (int j = 0; j < grid[0].size(); j++) {
        point_t point = copy_grid[row_index][j];
        copy_grid[row_index][j] = {point.x, point.t, f(point.x, point.t), true};
    }
    return copy_grid;
}

vector<vector<point_t>> fillGridRowTFixed(vector<vector<point_t>> grid, int column_index, double (*f)(double, double)) {
    vector<vector<point_t>> copy_grid(grid);
    for (int i = 0; i < grid.size(); i++) {
        point_t point = copy_grid[i][column_index];
        copy_grid[i][column_index] = {point.x, point.t, f(point.x, point.t), true};
    }
    return copy_grid;
}

vector<vector<point_t>> prepareGrid(
        interval_t x_domain,
        interval_t t_domain,
        int spatial_step_number,
        int time_step_number
) {
    vector<vector<point_t>> grid = getGrid(t_domain, x_domain, spatial_step_number, time_step_number);
    grid = fillGridRowXFixed(grid, 0, [](double x, double t) {
        return 0.0;
    });
    grid = fillGridRowXFixed(grid, grid.size() - 1, [](double x, double t) {
        return 0.0;
    });
    grid = fillGridRowTFixed(grid, 0, [](double x, double t) {
        return 0.1 * sin(M_PI * x);
    });
    grid = fillGridRowTFixed(grid, 1, [](double x, double t) {
        return 0.1 * sin(M_PI * x);
    });
    return grid;
}

//u_{i}^{j+1} = 2*u_i^j - u_i^{j-1} +
//      D_i * (h_t^2/h_x^2) * [u_{i+1}^j - 2*u_i^j + u_{i-1}^j]
double calculate_u_i_jpos(
        vector<vector<point_t>> grid,
        int i,
        int j,
        double spatial_step_number,
        double time_step_number,
        double(*D)(double)) {
    double u_i_j = grid[i][j].value;
    double u_i_jmin = grid[i][j - 1].value;
    double u_ipos_j = grid[i + 1][j].value;
    double u_imin_j = grid[i - 1][j].value;
    double ht = time_step_number;
    double hx = spatial_step_number;
    return 2 * u_i_j - u_i_jmin +
           D(i) * (ht * ht / hx * hx) * (u_ipos_j - 2 * u_i_j + u_imin_j);
}

double D(double x) {
    return exp(cos(x)) / 10;
}

vector<vector<point_t>> serial_approach(
        vector<vector<point_t>> grid_in,
        double spatial_step_number,
        double time_step_number
) {
    vector<vector<point_t>> grid(std::move(grid_in));
    for (int j = 0; j < grid[0].size(); j++) {
        if (j == 0 || j == grid[0].size() - 1) {
            continue;
        }
        for (int i = 0; i < grid.size(); i++) {
            if (i == 0 || i == grid.size() - 1) {
                continue;
            }
            point_t point = grid[i][j + 1];
            grid[i][j + 1] =
                    {point.x,
                     point.t,
                     calculate_u_i_jpos(
                             grid,
                             i,
                             j,
                             spatial_step_number,
                             time_step_number,
                             D),
                     true
                    };
        }
    }
    return grid;
}

void renderGrid(vector<vector<point_t>> grid) {
    FILE *pipe = popen("gnuplot -persist", "w");

    printf("1");
    if (pipe != NULL) {
        printf("2");
        fprintf(pipe, "set terminal png\n");
        fprintf(pipe, "set output \"2D map.png\"\n");
        fprintf(pipe, "set view map\n");
        fprintf(pipe, "set dgrid3d\n");
        fprintf(pipe, "set xrange[0:1]\n");
        fprintf(pipe, "set yrange[0:2]\n");
        fprintf(pipe, "set zrange[-0.25:1.0]\n");
        fprintf(pipe, "set samples 25\n");
        fprintf(pipe, "set isosamples 20\n");
        fprintf(pipe, "set pm3d map interpolate 0,0\n"); // interpolate 0, 0
        fprintf(pipe, "set palette rgbformulae 22,13,-31\n");
        fprintf(pipe, "set title \"2D map\"\n");
        fprintf(pipe, "splot \"-\" using 1:2:3 with pm3d\n");


        for (int i = 0; i < grid.size(); i++) {
            for (int j = 0; j < grid[0].size(); j++) {
                point_t point = grid[i][j];
                fprintf(pipe, "%f %f %f\n", point.x, point.t, point.value);
            }
        }

        printf("Hello");
        fprintf(pipe, "e\n");
        fflush(pipe);
        pclose(pipe);
    }
}

int main() {
    int spatial_step_number = 200;
    int time_step_number = 1000;
    cout << "hello";
    interval_t x_domain = {0, 1};
    interval_t t_domain = {0, 2};
    vector<vector<point_t>> serial_solution =
            serial_approach(
                    prepareGrid(x_domain, t_domain, spatial_step_number, time_step_number),
                    spatial_step_number,
                    time_step_number
            );
    renderGrid(
            prepareGrid(x_domain, t_domain, spatial_step_number, time_step_number));

    std::cout << "Hello, World!" << std::endl;
    return 0;
}





