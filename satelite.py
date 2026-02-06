# File: satellite.py
# Creation: Saturday January 23rd 2021
# Author: Arthur Dujardin
# ------
# Copyright (c) 2021 Arthur Dujardin


# Basic imports
from abc import ABC

# GNSS Tools
from gnsstime import to_gnsstime, gnsstime


class Satellite(ABC):
    def __init__(self, prn, toc):
        self._prn = prn
        self._toc = to_gnsstime(toc)

    @property
    def system(self):
        return None

    @property
    def prn(self):
        return self._prn

    @property
    def toc(self):
        return self._toc

    def position(self, *args, **kwargs):
        raise NotImplementedError("This method is currently not available. Make a PR if you wish to update gnsstools.")

    def __repr__(self):
        rep = f"{self.__class__.__name__}("
        rep += f"\n  system: {self.system}"
        rep += f"\n  prn: {self.prn}"
        rep == f"\n  toc: {self.toc}"
        for attr, value in self.__dict__.items():
            if attr[0] != "_":
                rep += f"\n  {attr}: {value:.6e}"
        rep += "\n)"
        return rep

# File: gps.py
# Creation: Sunday January 24th 2021
# Author: Arthur Dujardin
# ------
# Copyright (c) 2021 Arthur Dujardin


# Basic imports
import numpy as np
import pandas as pd

# GNSS Tools
from satellite import Satellite
from utils import camel2snake


class GPS(Satellite):

    def __init__(self, prn=None, toc=None,
                 sv_clock_bias=None, sv_clock_drift=None, sv_clock_drift_rate=None,
                 iode=None, crs=None, delta_n=None, m0=None,
                 cuc=None, e=None, cus=None, sqrt_a=None,
                 toe=None, cic=None, omega0=None, cis=None,
                 i0=None, crc=None, omega=None, omega_dot=None,
                 idot=None, l2_codes=None, gps_week=None, l2_pflag=None,
                 sv_acc=None, sv_health=None, tgd=None, iodc=None,
                 trans_time=None, fit_inter=None):
        super().__init__(prn=prn, toc=toc)
        # First row
        self.sv_clock_bias = sv_clock_bias
        self.sv_clock_drift = sv_clock_drift
        self.sv_clock_drift_rate = sv_clock_drift_rate
        # Second row
        self.iode = iode
        self.crs = crs
        self.delta_n = delta_n
        self.m0 = m0
        # Third row
        self.cuc = cuc
        self.e = e
        self.cus = cus
        self.sqrt_a = sqrt_a
        # Fourth row
        self.toe = toe
        self.cic = cic
        self.omega0 = omega0
        self.cis = cis
        # Fifth row
        self.i0 = i0
        self.crc = crc
        self.omega = omega
        self.omega_dot = omega_dot
        # Sixth row
        self.idot = idot
        self.l2_codes = l2_codes
        self.gps_week = gps_week
        self.l2_pflag = l2_pflag
        # Seventh row
        self.sv_acc = sv_acc
        self.sv_health = sv_health
        self.tgd = tgd
        self.iodc = iodc
        # Eighth row
        self.trans_time = trans_time
        self.fit_inter = fit_inter

    @property
    def system(self):
        return "G"

    def __repr__(self):
        rep = f"GPS("
        # First line
        rep += f"\n  system:               {self.system}"
        rep += f"\n  prn:                  {self.prn:d}"
        rep += f"\n  toc:                  {self.toc} [UTC] (Time Of Clock)"
        rep += f"\n  sv_clock_bias:       {self.sv_clock_bias: .6e} [s]"
        rep += f"\n  sv_clock_drift:      {self.sv_clock_drift: .6e} [s/s]"
        rep += f"\n  sv_clock_drift_rate: {self.sv_clock_drift_rate: .6e} [s/s2]"
        # Second line
        rep += f"\n  iode:                {self.iode: 13} (Issue Of Data, Ephemeris)"
        rep += f"\n  crs:                 {self.crs: .6e} [m]"
        rep += f"\n  delta_n:             {self.delta_n: .6e} [rad/s]"
        rep += f"\n  m0:                  {self.m0: .6e} [rad]"
        # Third line
        rep += f"\n  cuc:                 {self.cuc: .6e} [rad]"
        rep += f"\n  e:                   {self.e: .6e} (Eccentricity)"
        rep += f"\n  cus:                 {self.cus: .6e} [rad]"
        rep += f"\n  sqrt_a:              {self.sqrt_a: .6e} [sqrt(m)]"
        # Fourth line
        rep += f"\n  toe:                 {self.toe: .6e} [sec of GPS week] (Time Of Ephemeris)"
        rep += f"\n  cic:                 {self.cic: .6e} [rad]"
        rep += f"\n  omega0:              {self.omega0: .6e} [rad]"
        rep += f"\n  cis:                 {self.cis: .6e} [rad]"
        # Fifth line
        rep += f"\n  i0:                  {self.i0: .6e} [rad]"
        rep += f"\n  crc:                 {self.crc: .6e} [m]"
        rep += f"\n  omega:               {self.omega: .6e} [rad]"
        rep += f"\n  omega_dot:           {self.omega_dot: .6e} [rad/s]"
        # Sixth line
        rep += f"\n  idot:                {self.idot: .6e} [rad/s]"
        rep += f"\n  l2_codes:            {self.l2_codes: 13} (codes on L2 channel)"
        rep += f"\n  gps_week:            {self.gps_week: 13} (to go with TOE)"
        rep += f"\n  l2_pflag:            {self.l2_pflag: 13} (L2 P data flag)"
        # Seventh line
        rep += f"\n  sv_acc:              {self.sv_acc: .6e} [m]"
        rep += f"\n  sv_health:           {self.sv_health: .6e} (bits 17-22 w 3 sf 1)"
        rep += f"\n  tgd:                 {self.tgd: .6e} [s]"
        rep += f"\n  iodc:                {self.iodc: 13} (Issue Of Data, Clock)"
        # Eighth line
        rep += f"\n  trans_time:          {self.trans_time: .6e} [sec of GPS week] (e.g. derived from Z-count in Hand Over Word (HOW))"
        rep += f"\n  fit_inter:           {self.fit_inter: .6e} [hours] (Fit Interval in hours)"
        rep += f"\n)"
        return rep


# File: glonass.py
# Creation: Sunday January 24th 2021
# Author: Arthur Dujardin
# ------
# Copyright (c) 2021 Arthur Dujardin


from satellite import Satellite


class GLONASS(Satellite):
    
    def __init__(self, prn=None, toc=None,
                 sv_clock_bias=None, sv_rel_freq_bias=None, message_frame_time=None,
                 x=None, dx=None, dx2=None, health=None,
                 y=None, dy=None, dy2=None, freq_num=None,
                 z=None, dz=None, dz2=None, age_op_info=None):
        super().__init__(prn=prn, toc=toc)
        self.sv_clock_bias = sv_clock_bias
        self.sv_rel_freq_bias = sv_rel_freq_bias
        self.message_frame_time = message_frame_time
        # Second row
        self.x = x
        self.dx = dx
        self.dx2 = dx2
        self.health = health
        # Third row
        self.y = y
        self.dy = dy
        self.dy2 = dy2
        self.freq_num = freq_num
        # Fourth row
        self.z = z
        self.dz = dz
        self.dz2 = dz2
        self.age_op_info = age_op_info
        
    @property
    def system(self):
        return "R"

    def __repr__(self):
        rep = f"GLONASS("
        rep += f"\n  system:              {self.system}"
        rep += f"\n  prn:                 {self.prn:d}"
        rep += f"\n  toc:                 {self.toc} [UTC] (Time Of Clock)"
        rep += f"\n  sv_clock_bias:      {self.sv_clock_bias: .6e} [s] (-TauN)" 
        rep += f"\n  sv_rel_freq_bias:   {self.sv_rel_freq_bias: .6e} [s] (+GammaN)"
        rep += f"\n  message_frame_time: {self.message_frame_time: .6e} [s] (tk + nd * 86400 in seconds of the UTC week)"
        # Second line
        rep += f"\n  x:                  {self.x: .6e} [km] (satellite position X)"
        rep += f"\n  dx:                 {self.dx: .6e} [km/s] (velocity X dot)"
        rep += f"\n  dx2:                {self.dx2: .6e} [km/s2] (X acceleration)"
        rep += f"\n  health:             {self.health: 13} (0=healthy, 1=unhealthy)"
        # Third line
        rep += f"\n  y:                  {self.y: .6e} [km] (satellite position Y)"
        rep += f"\n  dy:                 {self.dy: .6e} [km/s] (velocity Y dot)"
        rep += f"\n  dy2:                {self.dy2: .6e} [km/s2] (Y acceleration)"
        rep += f"\n  freq_num:           {self.freq_num: 13} (frequency number)"
        # Fourth line
        rep += f"\n  z:                  {self.z: .6e} [km] (satellite position Z)"
        rep += f"\n  dz:                 {self.dz: .6e} [km/s] (velocity Z dot)"
        rep += f"\n  dz2:                {self.dz2: .6e} [km/s2] (Z acceleration)"
        rep += f"\n  age_op_info:        {self.age_op_info: 13} [days]"
        rep += f"\n)"
        return rep


# File: galileo.py
# Creation: Sunday January 24th 2021
# Author: Arthur Dujardin
# ------
# Copyright (c) 2021 Arthur Dujardin


from satellite import Satellite


class GALILEO(Satellite):

    def __init__(self, prn=None, toc=None,
                 sv_clock_bias=None, sv_clock_drift=None, sv_clock_drift_rate=None,
                 iod_nav=None, crs=None, delta_n=None, m0=None,
                 cuc=None, e=None, cus=None, sqrt_a=None,
                 toe=None, cic=None, omega0=None, cis=None,
                 i0=None, crc=None, omega=None, omega_dot=None,
                 idot=None, gps_week=None, gal_week=None,
                 sisa=None, sv_health=None, bgd_e5a=None, bgd_e5b=None,
                 trans_time=None):
        super().__init__(prn=prn, toc=toc)
        # First row
        self.sv_clock_bias = sv_clock_bias
        self.sv_clock_drift = sv_clock_drift
        self.sv_clock_drift_rate = sv_clock_drift_rate
        # Second row
        self.iod_nav = iod_nav
        self.crs = crs
        self.delta_n = delta_n
        self.m0 = m0
        # Third row
        self.cuc = cuc
        self.e = e
        self.cus = cus
        self.sqrt_a = sqrt_a
        # Fourth row
        self.toe = toe
        self.cic = cic
        self.omega0 = omega0
        self.cis = cis
        # Fifth row
        self.i0 = i0
        self.crc = crc
        self.omega = omega
        self.omega_dot = omega_dot
        # Sixth row
        self.idot = idot
        self.gps_week = gps_week
        self.gal_week = gal_week
        # Seventh row
        self.sisa = sisa
        self.sv_heath = sv_health
        self.bgd_e5a = bgd_e5a
        self.bgd_e5b = bgd_e5b
        # Eighth row
        self.trans_time = trans_time

    @property
    def system(self):
        return "E"

    def __repr__(self):
        rep = f"GALILEO("
        # First line
        rep += f"\n  system:               {self.system}"
        rep += f"\n  prn:                  {self.prn:d}"
        rep += f"\n  toc:                  {self.toc} [UTC] (Time Of Clock)"
        rep += f"\n  sv_clock_bias:       {self.sv_clock_bias: .6e} [s]"
        rep += f"\n  sv_clock_drift:      {self.sv_clock_drift: .6e} [s/s]"
        rep += f"\n  sv_clock_drift_rate: {self.sv_clock_drift_rate: .6e} [s/s2]"
        # Second line
        rep += f"\n  iod_nav:             {self.iod_nav: .6e} (Issue Of Data of the nav batch)"
        rep += f"\n  crs:                 {self.crs: .6e} [m]"
        rep += f"\n  delta_n:             {self.delta_n: .6e} [rad/s]"
        rep += f"\n  m0:                  {self.m0: .6e} [rad]"
        # Third line
        rep += f"\n  cuc:                 {self.cuc: .6e} [rad]"
        rep += f"\n  e:                   {self.e: .6e} (Eccentricity)"
        rep += f"\n  cus:                 {self.cus: .6e} [rad]"
        rep += f"\n  sqrt_a:              {self.sqrt_a: .6e} [sqrt(m)]"
        # Fourth line
        rep += f"\n  toe:                 {self.toe: .6e} [sec of GAL week] (Time Of Ephemeris)"
        rep += f"\n  cic:                 {self.cic: .6e} [rad]"
        rep += f"\n  omega0:              {self.omega0: .6e} [rad]"
        rep += f"\n  cis:                 {self.cis: .6e} [rad]"
        # Fifth line
        rep += f"\n  i0:                  {self.i0: .6e} [rad]"
        rep += f"\n  crc:                 {self.crc: .6e} [m]"
        rep += f"\n  omega:               {self.omega: .6e} [rad]"
        rep += f"\n  omega_dot:           {self.omega_dot: .6e} [rad/s]"
        # Sixth line
        rep += f"\n  idot:                {self.idot: .6e} [rad/s]"
        rep += f"\n  l2_codes:            {self.l2_codes: .6e} (codes on L2 channel)"
        rep += f"\n  gal_week:            {self.gal_week: .6e} (to go with TOE)"
        # Seventh line
        rep += f"\n  sisa:                {self.sisa: .6e} [m] (Signal in space accuracy)"
        rep += f"\n  sv_health:           {self.sv_health: .6e} (See Galileo ICD Section 5.1.9.3)"
        rep += f"\n  bgd_e5a:             {self.bgd_e5a: .6e} [s] (BGD E5a/E1)"
        rep += f"\n  bgd_e5b:             {self.bgd_e5b: .6e} [s] (BGD E5b/E1)"
        # Eighth line
        rep += f"\n  trans_time:          {self.trans_time: .6e} [sec of GAL week] (e.g. derived from WN and TOW of page type 1)"
        rep += f"\n)"
        return rep
