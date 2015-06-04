# -*- coding: utf-8 -*-
import matplotlib
import matplotlib.pyplot as plt


def strong_scaling():
  data = [
  # 64    128     256     512     1024    2048
  [2.81,  1.76,   1.08,   0.91,   1.01,    6.7],
  [6.15,  3.71,   2.01,   1.67,   1.33,    7.7],
  [13.00, 7.63,   3.76,   2.91,   2.58,    8.7],
  [       15.02,  7.78,   5.85,   4.53,    1.2],
  [               15.25,  10.95,  8.89,   10.6],
  [                       22.59,  17.77,   10.5],
  ]

  x_domain = [64, 128, 256, 512, 1024, 2048]
  legend_labels = ['Scale 28', 'Scale 29','Scale 30','Scale 31']

  f, ax = plt.subplots(1, 1)
  ax.plot(x_domain, data[0], '.-', linewidth=2.0, ms=10, label='Scale 26')
  ax.plot(x_domain, data[1], '.-', linewidth=2.0, ms=10, label='Scale 27')
  ax.plot(x_domain, data[2], '.-', linewidth=2.0, ms=10, label='Scale 28')
  ax.plot(x_domain[1:], data[3], '.-', linewidth=2.0, ms=10, label='Scale 29')
  ax.plot(x_domain[2:], data[4], '.-', linewidth=2.0, ms=10, label='Scale 30')
  ax.plot(x_domain[3:], data[5], '.-', linewidth=2.0, ms=10, label='Scale 31')
  ax.set_title('Stong Scaling')
  # ax.set_xticks([0,1,2,3,4])
  ax.set_ylim([0, 23])
  
  handles, labels = ax.get_legend_handles_labels()
  # sort both labels and handles by labels
  labels, handles = zip(*sorted(zip(labels, handles), key=lambda t: t[0], reverse=True))
  ax.legend(handles, labels, bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.)

  ax.set_xlim([58, 2300])
  ax.set_xscale('log')
  ax.set_ylabel('Time Per Stream (s)')
  ax.set_xlabel('# MPI Processes (log)')
  ax.set_xticks(x_domain)
  ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())

  plt.show()

def bipartite():
  biz_trans = {
  'amazon0302':[ 0.201, 0.133, 0.114, 0.107, 0.105, 0.104, 0.103, 0.103, 0.103, 0.102],
  'roadNet-CA.mtx':[ 0.186, 0.129, 0.112, 0.104, 0.098, 0.094, 0.091, 0.089, 0.088, 0.086] 
  }

  cite_collab = {
  'ca-AstroPh':[ 0.204, 0.145, 0.126, 0.115, 0.109, 0.104, 0.1, 0.095, 0.094, 0.093],
  'ca-GrQc':[ 0.136, 0.097, 0.088, 0.084, 0.082, 0.081, 0.081, 0.08, 0.08, 0.079],
  'ca-HepPh':[ 0.114, 0.075, 0.071, 0.07, 0.07, 0.07, 0.07, 0.07, 0.07, 0.07],
  'cit-Patents':[ 0.401, 0.343, 0.29, 0.256, 0.231, 0.241, 0.222, 0.216, 0.215, 0.209]}

  commune_social = {
  'email-EuAll':[0.27, 0.11, 0.101, 0.099, 0.097, 0.097, 0.097, 0.096, 0.081, 0.081],
  'soc-LiveJournal':[ 0.237, 0.156, 0.129, 0.112, 0.098, 0.085, 0.074, 0.064, 0.057, 0.051],
  'wiki-Talk':[ 0.412, 0.316, 0.241, 0.244, 0.189, 0.125, 0.088, 0.075, 0.071, 0.081] }

  as_web_related = {
  'Oregon-1':[ 0.266, 0.109, 0.095, 0.092, 0.091, 0.089, 0.088, 0.085, 0.083, 0.081],
  'web-Google':[ 0.184, 0.111, 0.095, 0.086, 0.091, 0.096, 0.096, 0.095, 0.092, 0.089],
  'p2p-Gnutella04':[ 0.412, 0.364, 0.306, 0.273, 0.254, 0.248, 0.245, 0.24, 0.24, 0.239],
  'as-Skitter':[ 0.171, 0.103, 0.093, 0.089, 0.087, 0.086, 0.085, 0.084, 0.084, 0.083]}



  params = {'legend.fontsize': 11,
          'legend.linewidth': 2}
  plt.rcParams.update(params)

  f, (ax1, ax2, ax3) = plt.subplots(3, 1) # ,sharey=True)

  lim = 5

  ax1.plot(cite_collab['ca-AstroPh'][0:lim], '.-')
  ax1.plot(cite_collab['ca-GrQc'][0:lim], '.-')
  ax1.plot(cite_collab['ca-HepPh'][0:lim], '.-')
  ax1.plot(cite_collab['cit-Patents'][0:lim], '.-')
  ax1.set_title('Citation & Collaboration')
  ax1.set_xticks([0,1,2,3,4])
  ax1.set_xticklabels(['1', '2', '3', '4', '5'])
  ax1.legend(['ca-AstroPh', 'ca-GrQc', 'ca-HepPh', 'cit-Patents'])
  
  ax2.plot(commune_social['email-EuAll'][0:lim], '.-')
  ax2.plot(commune_social['soc-LiveJournal'][0:lim], '.-')
  ax2.plot(commune_social['wiki-Talk'][0:lim], '.-')
  ax2.set_title('Communication & Social')
  ax2.set_ylabel('Fraction Edges Cut')
  ax2.set_xticks([0,1,2,3,4])
  ax2.set_xticklabels(['1', '2', '3', '4', '5'])
  ax2.legend(['email-EuAll', 'soc-LiveJournal', 'wiki-Talk'])

  ax3.plot(as_web_related['Oregon-1'][0:lim], '.-')
  ax3.plot(as_web_related['web-Google'][0:lim], '.-')
  ax3.plot(as_web_related['p2p-Gnutella04'][0:lim], '.-')
  ax3.plot(as_web_related['as-Skitter'][0:lim], '.-')
  ax3.legend(['Oregon-1', 'web-Google', 'p2p-Gnutella04', 'as-Skitter'])
  ax3.set_title('Autonomous Sys., Traffic, & Web')
  ax3.set_xticks([0,1,2,3,4])
  ax3.set_xticklabels(['1', '2', '3', '4', '5'])
  ax3.set_xlabel('Iterations')

  # plt.title('Fraction of Edges Cut (  ) for 2-Partitions')
  f.set_size_inches(7,9)
  f.savefig('real_k2_lambda.pdf', dpi=150)

  # plt.show()

def second_plot():





  biz_trans = {
  'amazon0302':[0.395,0.285,0.251,0.236,0.227,0.222,0.218,0.216,0.214,0.213],
  'roadNet-CA':[0.393,0.3,0.263,0.242,0.228,0.219,0.212,0.207,0.202,0.199]
  }

  cite_collab = {
  'ca-AstroPh':[0.353,0.323,0.287,0.226,0.219,0.216,0.209,0.208,0.208,0.207],
  'ca-GrQc':[0.269,0.189,0.175,0.17,0.167,0.166,0.165,0.165,0.165,0.164],
  'ca-HepPh':[0.232,0.195,0.185,0.182,0.181,0.181,0.18,0.18,0.18,0.179],
  'cit-Patents':[0.674,0.568,0.542,0.522,0.484,0.492,0.519,0.54,0.505,0.51]
  }

  commune_social = {
  'email-EuAll':[0.513,0.2,0.198,0.198,0.202,0.199,0.213,0.219,0.218,0.211],
  'soc-LiveJournal':[0.424,0.384,0.332,0.299,0.281,0.281,0.278,0.276,0.276,0.276],
  'wiki-Talk':[0.742,0.666,0.655,0.644,0.628,0.607,0.581,0.549,0.498,0.429]
  }

  as_web_related = {
  'Oregon-1':[0.415,0.302,0.291,0.285,0.283,0.281,0.281,0.281,0.28,0.28],
  'web-Google':[0.357,0.21,0.186,0.176,0.172,0.157,0.161,0.163,0.161,0.159],
  'p2p-Gnutella04':[0.815,0.763,0.763,0.759,0.754,0.748,0.748,0.737,0.732,0.732],
  'as-Skitter':[0.341,0.206,0.186,0.179,0.175,0.173,0.172,0.171,0.171,0.171]
  }



  params = {'legend.fontsize': 11,
          'legend.linewidth': 2}
  plt.rcParams.update(params)

  f, (ax1, ax2, ax3) = plt.subplots(3, 1) # ,sharey=True)

  lim = 5

  ax1.plot(cite_collab['ca-AstroPh'][0:lim], '.-')
  ax1.plot(cite_collab['ca-GrQc'][0:lim], '.-')
  ax1.plot(cite_collab['ca-HepPh'][0:lim], '.-')
  ax1.plot(cite_collab['cit-Patents'][0:lim], '.-')
  ax1.set_title('Citation & Collaboration')
  ax1.set_xticks([0,1,2,3,4])
  ax1.set_xticklabels(['1', '2', '3', '4', '5'])
  ax1.legend(['ca-AstroPh', 'ca-GrQc', 'ca-HepPh', 'cit-Patents'])
  
  ax2.plot(commune_social['email-EuAll'][0:lim], '.-')
  ax2.plot(commune_social['soc-LiveJournal'][0:lim], '.-')
  ax2.plot(commune_social['wiki-Talk'][0:lim], '.-')
  ax2.set_title('Communication & Social')
  ax2.set_ylabel('Fraction Edges Cut')
  ax2.set_xticks([0,1,2,3,4])
  ax2.set_xticklabels(['1', '2', '3', '4', '5'])
  ax2.legend(['email-EuAll', 'soc-LiveJournal', 'wiki-Talk'])

  ax3.plot(as_web_related['Oregon-1'][0:lim], '.-')
  ax3.plot(as_web_related['web-Google'][0:lim], '.-')
  ax3.plot(as_web_related['p2p-Gnutella04'][0:lim], '.-')
  ax3.plot(as_web_related['as-Skitter'][0:lim], '.-')
  ax3.legend(['Oregon-1', 'web-Google', 'p2p-Gnutella04', 'as-Skitter'])
  ax3.set_title('Autonomous Sys., Traffic, & Web')
  ax3.set_xticks([0,1,2,3,4])
  ax3.set_xticklabels(['1', '2', '3', '4', '5'])
  ax3.set_xlabel('Iterations')

  # plt.show()s
  # plt.title('Fraction of Edges Cut (  ) for 16-Partitions')
  f.set_size_inches(7,9)
  f.savefig('real_k16_lambda.pdf', dpi=150)

def scattar():
  terri = {
    'amazon0302':[0.00001797435054,  0.222],
    'amazon0601':[0.00002281642164,  0.198],
    'as-735':[0.0004445496569, 0.207],
    'as-Skitter':[0.00000771089446,  0.166],
    'ca-AstroPh':[0.001124215405,  0.232],
    'ca-CondMat':[0.0001493244869, 0.214],
    'ca-GrQc':[ 0.001054640264,  0.128],
    'ca-HepPh':[0.001643710433,  0.082],
    'ca-HepTh':[0.0005827346756, 0.202],
    'cit-HepPh':[ 0.0003532501881, 0.343],
    'cit-HepTh':[ 0.0004574940328, 0.36],
    'cit-Patents':[ 0.000001159316072, 0.402],
    'email-Enron':[ 0.0002730901121, 0.132],
    'email-EuAll':[ 0.000005971768011, 0.28],
    'Oregon-1':[0.0003545043941, 0.224],
    'Oregon-2':[0.0004696458004, 0.185],
    'p2p-Gnutella04':[0.0003379223282, 0.415],
    'roadNet-CA':[0.000001423902967, 0.186],
    'roadNet-PA':[0.000002591193426, 0.188],
    'roadNet-TX':[0.000001979545692, 0.172],
    'soc-Epinions1':[ 0.00008835527213,  0.173],
    'soc-LiveJournal1':[0.000002936037093, 0.234],
    'soc-sign-epinions':[ 0.00004841419648,  0.161],
    'soc-Slashdot0811':[0.0001513004377, 0.179],
    'soc-Slashdot0902':[0.0001404802977, 0.236],
    'web-BerkStan':[0.00001618731636, 0.187],
    'web-Google':[0.000006078583186, 0.189],
    'web-NotreDame':[ 0.00001411067604,  0.204],
    'web-Stanford':[0.00002909924913,  0.154],
    'wiki-Talk':[ 0.0000008758660987,  0.411],
    'wiki-Vote':[ 0.001506227269,  0.409]
  }
  data = [[],[]]
  for k in terri.keys():
    data[0].append(terri[k][0])
    data[1].append(terri[k][1])

  
  labels = terri.keys()
  plt.subplots_adjust(bottom = 0.1)
  f, ax1 = plt.subplots(1, 1) # ,sharey=True)
  plt.scatter(data[0], data[1], marker = 'o', color = '#4B6E9C')
  for label, x, y in zip(labels, data[0], data[1]):
      plt.annotate(
          label, 
          xy = (x, y), xytext = (0, -12), fontsize=10,
          textcoords = 'offset points', ha = 'center', va = 'bottom')

  # ax1.scatter(data[0], data[1])
  ax1.set_xlim([0.0000005, 0.0025])
  ax1.set_xscale('log')
  ax1.set_ylabel('Fraction Edges Cut (5 passes)')
  ax1.set_xlabel('Nonzero Density (log)')
  plt.show()

if __name__ == '__main__':
  # bipartite()
  # second_plot()
  # scattar()
  strong_scaling()