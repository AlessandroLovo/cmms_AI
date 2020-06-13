#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun  6 16:25:01 2020

@author: alessandro
"""

import os
import time
import pandas as pd
import numpy as np


pseudo_dir = '/home/alessandro/Quantum-Espresso/qe-6.5/pseudo/'
scratch_dir = '/home/alessandro/scratch'
pw_exe = '~/Quantum-Espresso/qe-6.5/bin/pw.x'

pseudo_name_list = [name for name in os.listdir(pseudo_dir) if name.startswith('Al')]

def make_bulk_scf_df_template():
    return pd.DataFrame(data=[],columns=['pseudo_name','a_0','ecutwfc','ecutrho_r','degauss','mixing_beta','conv_thr_o','n_k_points',
                        'cell_volume','pressure','total_energy','n_it','sim_time'])


def bulk_sim_scf(work_dir='.',name='al-bulk-test-scf', overwrite=True, verbose=False, save_df = None,
                 prefix='Al_bulk',
                 a_0=7.6,ecutwfc=60.0,ecutrho_r=4,degauss=0.05,
                 mixing_beta=0.7,conv_thr_o=-8,
                 pseudo_name='Al.pz-vbc.UPF',
                 fixed_atom='1 1 1',
                 n_k_points=6, k_points_shift=True):
    '''
    a_0 in Bohr
    energies in Ry
    press in kbar
    forces in hartree/bohr
    '''
    
    start_time = time.time()
    
    filename_in = work_dir.rstrip('/') + '/' + name + '.in'
    filename_out = work_dir.rstrip('/') + '/' + name + '.out'
    
    if os.path.exists(filename_in) and not overwrite:
        raise FileExistsError(filename_in + ' exists! To overwrite set overwrite=True')
    # write input for pw.x
    file = open(filename_in,'w')
    
    file.write('&CONTROL\n')
    file.write('  restart_mode=\'from_scratch\',\n')
    file.write('  prefix=\'%s\',\n' %prefix)
    file.write('  pseudo_dir = \'%s\',\n' %pseudo_dir)
    file.write('  outdir = \'%s\',\n' %scratch_dir)
    file.write('  calculation = \'scf\',\n')
    file.write('  tprnfor = .true.,\n')
    file.write('  tstress = .true.,\n')
    file.write('/\n')
    
    file.write('&SYSTEM\n')
    file.write('  ibrav = 2,\n')
    file.write('  celldm(1) = %f,\n' %a_0)
    file.write('  nat = 1,\n')
    file.write('  ntyp = 1,\n')
    file.write('  ecutwfc = %.1f,\n' %ecutwfc)
    file.write('  ecutrho = %.1f,\n' %(ecutwfc*ecutrho_r))
    file.write('  occupations = \'smearing\',\n')
    file.write('  smearing = \'mp\',\n')
    file.write('  degauss = %.4f,\n' %degauss)
    file.write('/\n')
    
    file.write('&ELECTRONS\n')
    file.write('  mixing_beta = %.2f,\n' %mixing_beta)
    file.write('  conv_thr = 1d%d,\n' %conv_thr_o)
    file.write('/\n')
    
    file.write('ATOMIC_SPECIES\n')
    file.write('  Al 26.981 %s\n' %pseudo_name)
    
    file.write('ATOMIC_POSITIONS (alat)\n')
    file.write('  Al 0. 0. 0. %s\n' %fixed_atom)
    
    file.write('K_POINTS automatic\n')
    if k_points_shift:
        file.write('  %d %d %d 1 1 1\n' %(n_k_points,n_k_points,n_k_points))
    else:
        file.write('  %d %d %d 0 0 0\n' %(n_k_points,n_k_points,n_k_points))
        
    file.close()
    
    #run pw.x
    cmd = '%s < %s > %s' %(pw_exe,filename_in,filename_out)
    if verbose:
        print(cmd)
    os.system('echo '+cmd)
    os.system(cmd)
    
    
    #collect cell volume
    os.system('grep \'unit-cell volume\' %s > tmp.txt' %filename_out)
    line = open('tmp.txt','r').readline()
    
    prefix,sep,suffix = line.partition('=')
    v_s,sep,suffix = suffix.lstrip(' ').partition(' ')
    v = float(v_s)
    
    # collect pressure
    os.system('grep \'P=\' %s > tmp.txt' %filename_out)
    line = open('tmp.txt','r').readline()
    prefix,sep,suffix = line.partition('=')
    p = float(suffix.lstrip(' '))

    #collect total energiy
    os.system('grep \'!    total energy\' %s > tmp.txt' %filename_out)
    line = open('tmp.txt','r').readline()
    prefix,sep,suffix = line.partition('=')
    e_s,sep,suffix = suffix.lstrip(' ').partition(' ')
    e = float(e_s)
    
    #collect number of iterations
    os.system('grep \'achieved\' %s > tmp.txt' %filename_out)
    line = open('tmp.txt','r').readline()
    try:
        n_it = int(line.split(' ')[-2])
    except(ValueError):
        print('Unable to compute number of iterations')
        n_it = 0
    
    os.remove('tmp.txt')
    t = time.time() - start_time
    
    if type(save_df) != type(None):
        save_df.loc[len(save_df)] = [pseudo_name,a_0,ecutwfc,ecutrho_r,degauss,mixing_beta,conv_thr_o,n_k_points,
                                        v,p,e,n_it,t]
    
    if verbose:
        print('Simulation took %.2f seconds' %t)
    return v, p, e, n_it, t



def bulk_sim_vc_relax(work_dir='.',name='al-bulk-test-vcr', overwrite=True, verbose=False,
                      prefix='Al_bulk',forc_conv_thr_o=-5,nstep=50,
                      a_0=7.6,ecutwfc=60.0,ecutrho_r=4,degauss=0.05,
                      mixing_beta=0.7,conv_thr_o=-8,
                      bfgs_ndim=1,upscale=100,pot_extrapolation='second_order',wfc_extrapolation='second_order',
                      press=0.0,press_conv_thr=0.5,
                      pseudo_name='Al.pz-vbc.UPF',
                      fixed_atom='1 1 1',
                      n_k_points=6, k_points_shift=True):
    '''
    a_0 in Bohr
    energies in Ry
    press in kbar
    forces in hartree/bohr
    '''
    
    start_time = time.time()
    
    filename_in = work_dir.rstrip('/') + '/' + name + '.in'
    filename_out = work_dir.rstrip('/') + '/' + name + '.out'
    
    if os.path.exists(filename_in) and not overwrite:
        raise FileExistsError(filename_in + ' exists! To overwrite set overwrite=True')
    # write input for pw.x
    file = open(filename_in,'w')
    
    file.write('&CONTROL\n')
    file.write('  restart_mode=\'from_scratch\',\n')
    file.write('  prefix=\'%s\',\n' %prefix)
    file.write('  pseudo_dir = \'%s\',\n' %pseudo_dir)
    file.write('  outdir = \'%s\',\n' %scratch_dir)
    file.write('  calculation = \'vc-relax\',\n')
    file.write('  forc_conv_thr = 1d%d,\n' %forc_conv_thr_o)
    file.write('  nstep = %d,\n' %nstep)
    file.write('/\n')
    
    file.write('&SYSTEM\n')
    file.write('  ibrav = 2,\n')
    file.write('  celldm(1) = %f,\n' %a_0)
    file.write('  nat = 1,\n')
    file.write('  ntyp = 1,\n')
    file.write('  ecutwfc = %.1f,\n' %ecutwfc)
    file.write('  ecutrho = %.1f,\n' %(ecutwfc*ecutrho_r))
    file.write('  occupations = \'smearing\',\n')
    file.write('  smearing = \'mp\',\n')
    file.write('  degauss = %.4f,\n' %degauss)
    file.write('/\n')
    
    file.write('&ELECTRONS\n')
    file.write('  mixing_beta = %.2f,\n' %mixing_beta)
    file.write('  conv_thr = 1d%d,\n' %conv_thr_o)
    file.write('/\n')
    
    file.write('&IONS\n')
    file.write('  ion_dynamics = \'bfgs\',\n')
    file.write('  bfgs_ndim = %d,\n' %bfgs_ndim)
    file.write('  upscale = %.1f,\n' %upscale)
    file.write('  pot_extrapolation = \'%s\',\n' %pot_extrapolation)
    file.write('  wfc_extrapolation = \'%s\',\n' %wfc_extrapolation)
    file.write('/\n')
    
    file.write('&CELL\n')
    file.write('  cell_dynamics = \'bfgs\',\n')
    file.write('  press = %.3f,\n' %press)
    file.write('  press_conv_thr = %.3f,\n' %press_conv_thr)
    file.write('/\n')
    
    file.write('ATOMIC_SPECIES\n')
    file.write('  Al 26.981 %s\n' %pseudo_name)
    
    file.write('ATOMIC_POSITIONS (alat)\n')
    file.write('  Al 0. 0. 0. %s\n' %fixed_atom)
    
    file.write('K_POINTS automatic\n')
    if k_points_shift:
        file.write('  %d %d %d 1 1 1\n' %(n_k_points,n_k_points,n_k_points))
    else:
        file.write('  %d %d %d 0 0 0\n' %(n_k_points,n_k_points,n_k_points))
        
    file.close()
    
    #run pw.x
    cmd = '%s < %s > %s' %(pw_exe,filename_in,filename_out)
    if verbose:
        print(cmd)
    os.system('echo '+cmd)
    os.system(cmd)
    
    
    #collect cell volumes
    os.system('grep \'unit-cell volume\' %s > tmp.txt' %filename_out)
    raw_v = [line for line in open('tmp.txt','r').readlines()]
    vs = []
    for line in raw_v[:-1]:
        prefix,sep,suffix = line.partition('=')
        v,sep,suffix = suffix.lstrip(' ').partition(' ')
        vs.append(float(v))
    
    #collect pressures
    os.system('grep \'P=\' %s > tmp.txt' %filename_out)
    raw_v = [line for line in open('tmp.txt','r').readlines()]
    ps = []
    for line in raw_v:
        prefix,sep,suffix = line.partition('=')
        ps.append(float(suffix.lstrip(' ')))
        
    #collect total energies
    os.system('grep \'!    total energy\' %s > tmp.txt' %filename_out)
    raw_v = [line for line in open('tmp.txt','r').readlines()]
    es = []
    for line in raw_v:
        prefix,sep,suffix = line.partition('=')
        e,sep,suffix = suffix.lstrip(' ').partition(' ')
        es.append(float(e))
    
    os.remove('tmp.txt')
    if verbose:
        print('Simulation took %.2f seconds' %(time.time() - start_time))
    return vs,ps,es



def surface_sim_relax(work_dir='.',name='al-surf-test', overwrite=True, verbose=False,
                      prefix='Al_surf',forc_conv_thr=0.001,
                      a_0=7.46834,vacuum_thickness=16.,n_layers=5,ecutwfc=60.0,ecutrho_r=4,degauss=0.02,
                      mixing_beta=0.7,conv_thr_o=-8,
                      bfgs_ndim=1,upscale=100,pot_extrapolation='second_order',wfc_extrapolation='second_order',
                      pseudo_name='Al.pz-vbc.UPF',
                      n_k_points_hor=8,n_k_points_ver=1,k_points_shift='1 1 0'):
    '''
    a_0 in Bohr
    energies in Ry
    press in kbar
    forces in hartree/bohr
    '''
    
    start_time = time.time()
    
    if n_layers %2 == 0:
        raise ValueError('Please use an odd number of layers.')
        
    aspect_ratio = (n_layers*0.5*a_0 + vacuum_thickness)/(a_0/np.sqrt(2))
    
    filename_in = work_dir.rstrip('/') + '/' + name + '.in'
    filename_out = work_dir.rstrip('/') + '/' + name + '.out'
    
    if os.path.exists(filename_in) and not overwrite:
        raise FileExistsError(filename_in + ' exists! To overwrite set overwrite=True')
    # write input for pw.x
    file = open(filename_in,'w')
    
    file.write('&CONTROL\n')
    file.write('  restart_mode=\'from_scratch\',\n')
    file.write('  prefix=\'%s\',\n' %prefix)
    file.write('  pseudo_dir = \'%s\',\n' %pseudo_dir)
    file.write('  outdir = \'%s\',\n' %scratch_dir)
    file.write('  calculation = \'relax\',\n')
    file.write('  forc_conv_thr = %f,\n' %forc_conv_thr)
    file.write('/\n')
    
    file.write('&SYSTEM\n')
    file.write('  nosym = .true.,\n')
    file.write('  ibrav = 6,\n')
    file.write('  celldm(1) = %f,\n' %(a_0/np.sqrt(2)))
    file.write('  celldm(3) = %f,\n' %aspect_ratio)
    file.write('  nat = %d,\n' %n_layers)
    file.write('  ntyp = 1,\n')
    file.write('  ecutwfc = %.1f,\n' %ecutwfc)
    file.write('  ecutrho = %.1f,\n' %(ecutwfc*ecutrho_r))
    file.write('  occupations = \'smearing\',\n')
    file.write('  smearing = \'mp\',\n')
    file.write('  degauss = %.4f,\n' %degauss)
    file.write('/\n')
    
    file.write('&ELECTRONS\n')
    file.write('  mixing_mode = \'local-TF\',\n')
    file.write('  mixing_beta = %.2f,\n' %mixing_beta)
    file.write('  conv_thr = 1d%d,\n' %conv_thr_o)
    file.write('/\n')
    
    file.write('&IONS\n')
    file.write('  ion_dynamics = \'bfgs\',\n')
    file.write('  bfgs_ndim = %d,\n' %bfgs_ndim)
    file.write('  upscale = %.1f,\n' %upscale)
    file.write('  pot_extrapolation = \'%s\',\n' %pot_extrapolation)
    file.write('  wfc_extrapolation = \'%s\',\n' %wfc_extrapolation)
    file.write('/\n')
    
    
    file.write('ATOMIC_SPECIES\n')
    file.write('  Al 26.981 %s\n' %pseudo_name)
    
    file.write('ATOMIC_POSITIONS (alat)\n')
    for i in range(n_layers):
        if i == n_layers // 2:
            file.write('  Al 0. 0. %f 0 0 0\n' %(i/np.sqrt(2)))
        elif i % 2 == 0:
            file.write('  Al 0. 0. %f\n' %(i/np.sqrt(2)))
        else:
            file.write('  Al 0.5 0.5 %f\n' %(i/np.sqrt(2)))
    
    file.write('K_POINTS automatic\n')
    file.write('  %d %d %d %s\n' %(n_k_points_hor,n_k_points_hor,n_k_points_ver,k_points_shift))
        
    file.close()
    
    #run pw.x
    cmd = '%s < %s > %s' %(pw_exe,filename_in,filename_out)
    if verbose:
        print(cmd)
    os.system('echo '+cmd)
    os.system(cmd)
    
    #collect atomic positions
    os.system('grep \'Al    \' %s > tmp.txt' %filename_out)
    raw_v = [line for line in open('tmp.txt','r').readlines()]
    atomic_coords = []
    for i,line in enumerate(raw_v[1:]):
        if i % n_layers == 0:
            atomic_coords.append([])
        epurated = line.rstrip('0 \n').lstrip('Al ')
        coords = [float(s) for s in epurated.split(' ') if len(s) > 0]
        atomic_coords[-1].append(coords)
    
    #collect total energies
    os.system('grep \'!    total energy\' %s > tmp.txt' %filename_out)
    raw_v = [line for line in open('tmp.txt','r').readlines()]
    es = []
    for line in raw_v:
        prefix,sep,suffix = line.partition('=')
        e,sep,suffix = suffix.lstrip(' ').partition(' ')
        es.append(float(e))
    
    os.remove('tmp.txt')
    t = time.time() - start_time
    
    if verbose:
        print('Simulation took %.2f seconds' %t)
        
    return np.array(es), np.array(atomic_coords), t