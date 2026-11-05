"""mef_masonry — elementos finitos térmicos e termomecânicos para alvenaria de
blocos (BCEPQ, cerâmicos estruturais e de vedação, concreto), com furos redondos,
células retangulares ou furos horizontais, graute, junta seca ou de argamassa,
em Python puro
(NumPy, SciPy, Matplotlib).

Módulos:
  parametros  propriedades, geometria e protocolo de ensaio (com as fontes)
  malha       seção 2D com furos (malha em O), células retangulares ou seção vertical
              de furos horizontais e prisma 3D por extrusão (W6, H8, W15), com
              graute, argamassa e furos horizontais fiada a fiada
  elementos   matrizes de elemento rápidas (triângulo T3 e cunha W6)
  elementos3d elementos isoparamétricos genéricos W6, H8 e W15
  naolinear   dano de Mazars no bloco, na argamassa e no graute; juntas secas com folgas e atrito;
              fluência transitória (Anderberg e Thelandersson) no aquecimento sob carga
  parede      parede sob carga e incêndio: vigas corrotacionais com seção de fibras e fluência transitória
  dados_parede curvas do prisma BCEPQ (30 a 1000 °C) usadas como lei da parede
  sensibilidade análise de sensibilidade global da parede no ISO 834 (método de Morris)
  biblioteca  blocos, materiais, argamassas e grautes com fonte (preenchem o projeto)
  projeto     projeto do programa (esquema, validação, execução)
  visual      exportação para o pós-processador WebGL
  termico     condução transiente com convecção e radiação
  mecanico    termoelasticidade 3D e estado plano generalizado 2D
  curvas      curvas de forno e ISO 834
  estudos     calibração das juntas, ensaio da UFV em 3D, estudos da seção
  saida       VTK (ParaView), CSV e figuras
"""
__version__ = "1.0.0"
