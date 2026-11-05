@echo off
setlocal
cd /d "%~dp0"
title MEF_Masonry
echo ==================================================
echo   MEF_Masonry - alvenaria de blocos em temperatura
echo ==================================================
echo.
echo [1/3] Procurando o Python...
set "PY="
where py >nul 2>&1
if not errorlevel 1 (
  py --version
  if not errorlevel 1 set "PY=py"
)
if not defined PY (
  where python >nul 2>&1
  if not errorlevel 1 (
    python --version
    if not errorlevel 1 set "PY=python"
  )
)
if not defined PY goto sem_python
echo.
echo [2/3] Conferindo as bibliotecas numpy, scipy e matplotlib...
%PY% -c "import numpy, scipy, matplotlib" >nul 2>&1
if not errorlevel 1 goto iniciar
echo Instalando as bibliotecas. Isso so acontece na primeira vez e leva alguns minutos...
%PY% -m pip install -r requirements.txt
if errorlevel 1 (
  echo Tentando de novo na pasta do usuario...
  %PY% -m ensurepip --upgrade >nul 2>&1
  %PY% -m pip install --user -r requirements.txt
)
%PY% -c "import numpy, scipy, matplotlib"
if errorlevel 1 goto erro_pip
:iniciar
echo.
echo [3/3] Abrindo o programa no navegador.
echo Deixe esta janela aberta enquanto usa o programa. Para encerrar, feche esta janela.
echo.
%PY% app.py %*
if errorlevel 1 goto erro_app
goto fim

:sem_python
echo.
echo Nao encontrei o Python neste computador.
echo  1. Se voce baixou o python-manager .msix, de dois cliques nele e clique em Instalar.
echo     Depois abra o Prompt de Comando e rode:  py install 3.13
echo  2. Ou instale pelo site python.org, marcando "Add python.exe to PATH".
echo Depois disso, feche esta janela e de dois cliques de novo neste arquivo.
echo.
pause
goto fim

:erro_pip
echo.
echo Nao foi possivel instalar as bibliotecas numpy, scipy e matplotlib.
echo Confira a conexao com a internet e tente de novo. Se continuar,
echo copie as mensagens acima e envie para o suporte.
echo.
pause
goto fim

:erro_app
echo.
echo O programa parou com erro. Copie as mensagens acima e envie para o suporte.
echo.
pause

:fim
endlocal
