/*
 * swmp.h - C interface for SWMP (Surface Wave Multipathing)
 *
 * This file is part of swmp.
 * Copyright (C) 2023 CSIRO
 *
 * swmp is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 */

#ifndef SWMP_H
#define SWMP_H

#ifdef __cplusplus
extern "C" {
#endif

/*=============================================================================
 * Error Codes
 *===========================================================================*/

#define SWMP_SUCCESS 0
#define SWMP_ERROR_FILE -1
#define SWMP_ERROR_MEMORY -2
#define SWMP_ERROR_INVALID_PARAM -3

/*=============================================================================
 * Error Handling
 *===========================================================================*/

/**
 * Get last error information
 * @param code Output error code
 * @param message Output error message buffer
 * @param length Length of message buffer
 * @return 0 on success
 */
int get_last_error(int* code, char* message, int length);

/**
 * Clear error state
 * @return 0 on success
 */
int clear_error(void);

/*=============================================================================
 * Model Generation (modgen module)
 *===========================================================================*/

/**
 * Read model generation configuration
 * @param filename Configuration file path
 * @param length Length of filename string
 * @return 0 on success, negative on error
 */
int gen2dv_read_conf(const char* filename, int length);

/**
 * Generate 2D velocity model
 * @return 0 on success, negative on error
 */
int gen2dv_run(void);

/*=============================================================================
 * Observation Generation (pred2obs module)
 *===========================================================================*/

/**
 * Read observation generation configuration
 * @param filename Configuration file path
 * @param length Length of filename string
 * @return 0 on success, negative on error
 */
int creobs_read_conf(const char* filename, int length);

/**
 * Create observations from predictions (add noise)
 * @return 0 on success, negative on error
 */
int creobs_run(void);

/*=============================================================================
 * Ray Tracing - Configuration (rat module)
 *===========================================================================*/

/**
 * Read ray tracing configuration
 * Loads velocity model, sources, receivers, and parameters
 * @param filename Configuration file path
 * @param length Length of filename string
 * @return 0 on success, negative on error
 */
int rat_read_conf(const char* filename, int length);

/**
 * Get number of sources
 * @param n Output number of sources
 * @return 0 on success
 */
int get_number_of_sources(int* n);

/**
 * Enable or disable file output
 * @param enable 0=no file writes, 1=write files (backward compat)
 * @return 0 on success
 */
int set_enable_file_output(int enable);

/*=============================================================================
 * Ray Tracing - Execution
 *===========================================================================*/

/**
 * Execute wavefront tracking for a single source (parallel-safe)
 *
 * This function is designed for parallel execution. Each worker process
 * can call this independently with a different source_id. Results remain
 * in module-level memory for retrieval via get_*_from_memory() functions.
 *
 * @param source_id Source index (1-based)
 * @return 0 on success, negative on error
 */
int forward_single_source(int source_id);

/**
 * Execute wavefront tracking for all sources (sequential)
 *
 * This provides backward compatibility. For parallel execution,
 * use forward_single_source() from multiple processes.
 *
 * @return 0 on success, negative on error
 */
int forward(void);

/**
 * Resample velocity model to finer grid
 * @param nry Y resampling factor
 * @param nrx X resampling factor
 * @return 0 on success
 */
int resample_model(int nry, int nrx);

/*=============================================================================
 * Ray Tracing - Parameter Getters/Setters
 *===========================================================================*/

int get_dt(float* dt);
int set_dt(float dt);

int get_maxit(int* maxit);
int set_maxit(int maxit);

int get_nsenod(int* nsenod);
int set_nsenod(int nsenod);

int get_mode(int* mode);
int set_mode(int mode);

int get_solver(int* solver);
int set_solver(int solver);

int get_interp(int* interp);
int set_interp(int interp);

int get_mar(int* mar);
int set_mar(int mar);

int get_velint(int* velint);
int set_velint(int velint);

int get_rearth(float* rearth);
int set_rearth(float rearth);

int get_recmode(int* recmode);
int set_recmode(int recmode);

int set_do_rays(int flag);
int set_do_frechet(int flag);

/*=============================================================================
 * Ray Tracing - Model Data Access
 *===========================================================================*/

/**
 * Get velocity model metadata
 * @param x0 Output X origin
 * @param y0 Output Y origin
 * @param nx Output number of X grid points
 * @param ny Output number of Y grid points
 * @param dx Output X grid spacing
 * @param dy Output Y grid spacing
 * @param cn Output number of cushion nodes
 * @return 0 on success
 */
int get_model_meta_data(float* x0, float* y0, int* nx, int* ny,
                        float* dx, float* dy, int* cn);

/**
 * Get velocity model values (flattened)
 * @param val Output array of size (nx+cn*2)*(ny+cn*2)
 * @param nx Number of X grid points
 * @param ny Number of Y grid points
 * @param cn Number of cushion nodes
 * @return 0 on success
 */
int get_model_vector(float* val, int nx, int ny, int cn);

/**
 * Set velocity model values (flattened)
 * @param val Input array of size (nx+cn*2)*(ny+cn*2)
 * @param nx Number of X grid points
 * @param ny Number of Y grid points
 * @param cn Number of cushion nodes
 * @return 0 on success
 */
int set_model_vector(const float* val, int nx, int ny, int cn);

/**
 * Get resampled model metadata
 */
int get_resampled_model_meta_data(float* x0, float* y0, int* nx, int* ny,
                                   float* dx, float* dy, int* cn);

/**
 * Get resampled model values (flattened)
 */
int get_resampled_model_vector(float* val, int nx, int ny, int cn);

/*=============================================================================
 * Ray Tracing - Arrival Data (File-Based)
 *===========================================================================*/

/**
 * Read observed arrivals from file
 * @param filename File path
 * @param length Length of filename string
 * @return 0 on success
 */
int read_observations(const char* filename, int length);

/**
 * Read predicted arrivals from file
 * @param filename File path
 * @param length Length of filename string
 * @return 0 on success
 */
int read_predictions(const char* filename, int length);

/**
 * Get number of observations
 * @param n Output count
 * @return 0 on success
 */
int get_number_of_observations(int* n);

/**
 * Get number of predictions
 * @param n Output count
 * @return 0 on success
 */
int get_number_of_predictions(int* n);

/**
 * Get observations (Nx6 array: src, rec, arr_num, time, azi, unc)
 * @param val Output array
 * @param n Number of observations
 * @return 0 on success
 */
int get_observations(float* val, int n);

/**
 * Get predictions (Nx5 array: src, rec, arr_num, time, azi)
 * @param val Output array
 * @param n Number of predictions
 * @return 0 on success
 */
int get_predictions(float* val, int n);

/*=============================================================================
 * Ray Tracing - Arrival Data (In-Memory - NEW)
 *===========================================================================*/

/**
 * Get number of arrivals from current memory state
 *
 * After calling forward_single_source(), this returns the number of
 * arrivals computed for that source. Handles variable arrivals per receiver.
 *
 * @param n Output number of arrivals
 * @return 0 on success
 */
int get_number_of_arrivals_in_memory(int* n);

/**
 * Extract arrivals from memory as flat arrays
 *
 * Gets arrival data from the current module state (after forward_single_source).
 * Arrays are flattened to handle variable arrival counts per receiver.
 *
 * @param receivers Output receiver IDs (1-based)
 * @param arrival_nums Output arrival numbers (1-based)
 * @param times Output travel times
 * @param azis Output azimuths
 * @param spfs Output spreading factors
 * @param n Number of arrivals
 * @return 0 on success
 */
int get_arrivals_from_memory(int* receivers, int* arrival_nums,
                              float* times, float* azis, float* spfs, int n);

/*=============================================================================
 * Ray Tracing - Raypath Data (In-Memory - NEW)
 *===========================================================================*/

/**
 * Count raypaths in memory
 * @param n_paths Output number of raypaths
 * @param total_points Output total number of points across all raypaths
 * @return 0 on success
 */
int get_raypath_count_in_memory(int* n_paths, int* total_points);

/**
 * Get raypath metadata
 * @param receivers Output receiver IDs
 * @param arrival_nums Output arrival numbers
 * @param npts Output number of points per raypath
 * @param n Number of raypaths
 * @return 0 on success
 */
int get_raypath_metadata_from_memory(int* receivers, int* arrival_nums,
                                      int* npts, int n);

/**
 * Get raypath positions (flattened)
 * @param positions Output array of shape (total_points, 2)
 * @param total_points Total number of points
 * @return 0 on success
 */
int get_raypath_positions_from_memory(float* positions, int total_points);

/*=============================================================================
 * Ray Tracing - Jacobian (Frechet Matrix)
 *===========================================================================*/

/**
 * Read Jacobian from files
 * @return 0 on success
 */
int read_jacobian(void);

/**
 * Get sparse Jacobian dimensions
 * @param nr Output number of rows
 * @param nc Output number of columns
 * @param nnz Output number of non-zeros
 * @return 0 on success
 */
int get_sparse_jacobian_size(int* nr, int* nc, int* nnz);

/**
 * Get sparse Jacobian in CRS format
 * @param jrow Row indices
 * @param jcol Column indices
 * @param jval Values
 * @param n Number of non-zeros
 * @return 0 on success
 */
int get_sparse_jacobian(float* jrow, float* jcol, float* jval, int n);

/*=============================================================================
 * Ray Tracing - File Paths (Backward Compatibility)
 *===========================================================================*/

int get_arrival_prediction_filepath(char* fn, int length);
int get_raypath_prediction_filepath(char* fn, int length);
int get_wavefront_prediction_filepath(char* fn, int length);

#ifdef __cplusplus
}
#endif

#endif /* SWMP_H */
