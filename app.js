/* LAMPrime — desenho de primers LAMP no navegador (algoritmo real)
   Tm nearest-neighbor (SantaLucia 1998) + correção de sal e equivalente
   monovalente de Mg2+ (von Ahsen 2001); GC%; complexidade linguística (LC%);
   estabilidade ΔG das extremidades (critério Eiken ≤ -4 kcal/mol);
   enumeração das regiões F3/F2/F1/B1/B2/B3 com restrições de distância
   (PrimerExplorer V5 / NEB) e montagem de FIP/BIP/LF/LB.
   Estruturas secundárias (hairpin/homo/heterodímero) @63 °C re-ranqueiam os finalistas.
   i18n PT/EN embutido (botão de idioma), sem dependências externas.
   Refs: Notomi 2000; Nagamine 2002; Tomita 2008; SantaLucia 1998; von Ahsen 2001. */
(function() {
  const qs = (sel, el=document) => el.querySelector(sel);
  const qsa = (sel, el=document) => Array.from(el.querySelectorAll(sel));

  // ===================================================================
  // i18n — dicionário PT/EN. Strings: aplicadas via [data-i18n*]. Funções:
  // usadas no render dinâmico dos resultados (interpolação).
  // ===================================================================
  const I18N = {
    pt: {
      'tab.input':'Entrada', 'tab.params':'Parâmetros', 'tab.results':'Resultados',
      'seq.h2':'Sequência alvo',
      'seq.help':'Cole a sequência (A, T, G, C) ou informe um identificador (ex.: GenBank).',
      'id.label':'Identificador (opcional)', 'id.ph':'Ex.: NC_045512, MN908947, etc.',
      'seq.label':'Sequência', 'seq.ph':'Cole aqui a sequência de DNA...',
      'btn.fasta':'Carregar FASTA', 'btn.clear':'Limpar', 'btn.run':'Gerar primers',
      'refs.h3':'Referências e Publicações',
      'refs.lead':'Referências-chave implementadas no algoritmo do LAMPrime:',
      'refs.tools':'Regras de geometria dos primers: Eiken PrimerExplorer V5 (https://primerexplorer.jp/) e New England Biolabs — LAMP Primer Design (https://lamp.neb.com/).',
      'refs.github':'Código-fonte aberto (MIT): <a href="https://github.com/luizdione/LAMPrime_browser-based" target="_blank" rel="noopener" style="color:var(--accent);">github.com/luizdione/LAMPrime_browser-based</a>',
      'params.basic.h2':'Configurações básicas',
      'tm.inner.label':'Tm alvo (F1c/B1c/LF/LB)', 'tm.outer.label':'Tm alvo (F2/B2/F3/B3)',
      'size.f3b3.label':'Tamanho dos primers (F3/B3)', 'size.lflb.label':'Tamanho dos primers (LF/LB)',
      'loopseq.label':'Linker do loop interno (FIP/BIP) — 5´TTTTTT 3´',
      'fip.label':'Tamanho total do FIP (F1c-loop-F2)', 'bip.label':'Tamanho total do BIP (B1c-loop-B2)',
      'gc.label':'GC esperado',
      'gc.opt.4060':'40–60%', 'gc.opt.4065':'40–65% (teto canônico Eiken)', 'gc.opt.5060':'50–60% (ótimo)',
      'gc.opt.3565':'35–65%', 'gc.opt.4555':'45–55%',
      'loops.label':'Loop primers',
      'loops.opt.auto':'Automático', 'loops.opt.on':'Forçar inclusão', 'loops.opt.off':'Não incluir',
      'lc.label':'Complexidade linguística mínima (LC%)',
      'lc.help':'A complexidade linguística (LC%) mede a “riqueza de vocabulário” de uma sequência genética, a partir do número de combinações distintas de nucleotídeos (k-mers) observadas em relação ao total teoricamente possível, expresso em porcentagem (100% = máximo). Intervalo permitido: 70–90%. Valor padrão: 75%.',
      'react.h2':'Condições de reação',
      'pen.h2':'Penalidades e ranking',
      'pen.f1cb1c.label':'Limite F1c/B1c', 'pen.f2b2.label':'Limite F2/B2', 'pen.f3b3.label':'Limite F3/B3',
      'maxsets.label':'Número máximo de conjuntos a testar',
      'pen.help':'As penalidades servem apenas para ranquear os conjuntos de primers. Nenhum conjunto é reprovado por penalidade: os limites são referência de pontuação.',
      'geo.h2':'Geometria dos alvos',
      'geo.f2b2.label':'Amplicon F2–B2 (5′F2→5′B2)',
      'geo.opt.120180':'120–180 bp', 'geo.opt.120160':'120–160 bp (PrimerExplorer/NEB)', 'geo.opt.120200':'120–200 bp',
      'geo.f3f2.label':'Gap F3–F2 (3′F3→5′F2)',
      'geo.opt.060':'0–60 bp', 'geo.opt.020':'0–20 bp (Eiken loopamp)', 'geo.opt.040':'0–40 bp',
      'geo.f1b1.label':'Distância F1–B1', 'geo.f2f1.label':'Distância F2–F1', 'geo.b2b1.label':'Distância B2–B1',
      'geo.help':'Geometria conforme o padrão Eiken PrimerExplorer V5 / NEB: <b>amplicon F2–B2 (5′F2→5′B2) = 120–180 bp</b> (região amplificada); <b>gap F3–F2 (3′F3→5′F2) = 0–60 bp</b>; loops F2–F1 e B2–B1 = 40–60 bp; centro F1–B1 = 0–60 bp.',
      'struct.h2':'Estruturas secundárias (hibridação entre oligos)',
      'struct.treact':'Temperatura de reação', 'struct.hpth':'Limiar ΔG — hairpin',
      'struct.dimerth':'Limiar ΔG — dímeros (homo/hetero)', 'struct.enddg':'Estabilidade de extremidade (ΔG ≤)',
      'struct.help':'O LAMPrime avalia, à temperatura de reação, as estruturas que prejudicam a LAMP: hairpins (intramoleculares) e dímeros — hibridação entre oligos, homo e heterodímeros — pelo modelo nearest-neighbor (SantaLucia 1998). Conjuntos com qualquer estrutura mais estável (ΔG mais negativo) que o limiar são penalizados e descem no ranking; nenhum é reprovado. A estabilidade das extremidades de iniciação segue o critério Eiken (ΔG ≤ −4 kcal/mol por padrão).',
      'btn.restore':'Restaurar padrão', 'btn.apply':'Aplicar parâmetros',
      'img.caption':'Esquema ilustrativo da disposição dos primers no método LAMP.',
      'res.h2':'Resultados',
      'res.empty':'Ainda não há resultados. Gere primers para visualizá-los aqui.',
      'footer.html':'<strong>LAMPrime</strong> — desenho de primers LAMP no navegador (Tm nearest-neighbor SantaLucia 1998; geometria PrimerExplorer V5/NEB). Código aberto (MIT): <a href="https://github.com/luizdione/LAMPrime_browser-based" target="_blank" rel="noopener" style="color:var(--accent);">github.com/luizdione/LAMPrime_browser-based</a>.',
      'footer.lab':'Laboratório de Genética Molecular — IFRJ. De-Melo, Luiz Dione Barbosa — Copyright 2025–2026.',
      // dinâmicas
      'btn.applying':'Aplicado!', 'btn.generating':'Gerando...',
      'alert.lc':'Informe um valor de LC% entre 70 e 90.',
      'alert.shortSeq':'Sequência curta. Insira ao menos ~150–250 nt para um design LAMP confiável.',
      'alert.noSets':'Nenhum conjunto LAMP satisfez as restrições (Tm/distâncias). Tente uma sequência mais longa (~200–300 nt), alargue as faixas de Tm ou o GC.',
      'alert.designErr':'Erro no desenho: ',
      'res.summary':'Resumo do desenho', 'res.target':'Alvo', 'res.sets':'Conjuntos',
      'res.targetTm':'Tm alvo', 'res.salts':'Sais',
      'res.method':'Tm por nearest-neighbor (SantaLucia 1998) com correção de sal + equivalente de Mg²⁺ (von Ahsen 2001). A penalidade apenas ranqueia; avisos não reprovam o conjunto.',
      'res.outer':'Primers externos',
      'res.inner':'Primers internos  (FIP = F1c–loop–F2 · BIP = B1c–loop–B2)',
      'res.loopHdr':'Loop primers (aceleram a reação)',
      'res.positions':'Posições (1-based)', 'res.distances':'Distâncias',
      'res.struct': t => `Estrutura 2ª (@${t} °C)`, 'res.worstHp':'pior hairpin ΔG', 'res.worstDi':'pior dímero ΔG',
      'res.allOk':'✓ todos os primers dentro das faixas de Tm/GC/LC', 'res.copy':'copiar',
      'res.set':'Conjunto', 'res.penalty':'penalidade', 'res.amplicon':'amplicon',
      'warn.hairpin':'hairpin', 'warn.selfdimer':'self-dímero', 'warn.dimer':'dímero',
      'notes.geometry': g => `Geometria (padrão Eiken/NEB): amplicon F2–B2 ${g.f2b2[0]}–${g.f2b2[1]} nt · gap F3–F2 ${g.gap[0]}–${g.gap[1]} nt · loop F2–F1 ${g.loopF2F1[0]}–${g.loopF2F1[1]} nt · centro F1–B1 ${g.f1b1[0]}–${g.f1b1[1]} nt.`,
      'res.setsOf': n => `(de ${n} candidatos válidos)`,
      // especificidade
      'tab.spec':'Especificidade',
      'spec.h2':'Especificidade (triagem in silico)',
      'spec.help':'Etapa opcional: tria os primers do conjunto melhor ranqueado (Conjunto #1) contra sequências de fundo, sinalizando primers que poderiam anelar e iniciar em organismos próximos. Usa semente 3′ exata + extensão (k-mer/alinhamento), não ML.',
      'spec.acc.label':'Acessos NCBI (opcional)',
      'spec.acc.ph':'Ex.: NC_043815.1, AY648786.1 (vírgula/espaço)',
      'spec.examples':'Exemplos de reatividade cruzada:',
      'spec.examplesHelp':'Preenche o acesso NCBI do organismo congenérico; gere primers, busque o FASTA e clique em Triar.',
      'btn.fetchncbi':'Buscar FASTA no NCBI',
      'spec.bg.label':'Sequências de fundo (FASTA)',
      'spec.bg.ph':'Cole um ou mais FASTA de fundo (genomas/genes de organismos próximos)...',
      'spec.seed.label':'Âncora 3′ exata (semente)',
      'spec.maxmm.label':'Máx. de mismatches (primer inteiro)',
      'btn.screen':'Triar especificidade',
      'spec.privacy':'Privacidade: a triagem roda no seu navegador e os primers não saem da página. A rede só é acessada se você usar “Buscar FASTA no NCBI” (são enviados apenas os identificadores públicos).',
      'spec.res.h2':'Resultado da triagem',
      'spec.empty':'Sem triagem ainda. Gere primers (aba Entrada) e informe sequências de fundo.',
      'btn.screening':'Triando...', 'btn.fetching':'Buscando...',
      'alert.specNoDesign':'Gere primers primeiro (aba Entrada → Gerar primers). A triagem usa o Conjunto #1.',
      'alert.specNoBg':'Informe ao menos uma sequência de fundo (cole o FASTA ou busque por acesso no NCBI).',
      'alert.specNoAcc':'Digite ao menos um número de acesso (ex.: NC_043815.1).',
      'alert.specFetchErr':'Falha ao buscar no NCBI (CORS ou rede). Baixe o FASTA manualmente no NCBI e cole no campo. Detalhe: ',
      'spec.bgHdr': (name, len) => `Fundo: ${name} (${len} nt)`,
      'spec.setEval': r => `Conjunto avaliado: #${r}`,
      'spec.summary': (hit, tot) => `${hit} de ${tot} primers com possível anelamento 3′-ancorado`,
      'spec.hit': (n, mm, pos, strand) => `⚠ ${n}: ${mm} mismatch(es), pos ${pos} (fita ${strand})`,
      'spec.clean': n => `✓ ${n}: sem anelamento 3′-ancorado`,
      'spec.multi':'Atenção: vários primers anelando ao mesmo fundo podem indicar risco de amplificação cruzada.',
      'spec.note':'Um hit é o melhor alinhamento sem gaps com semente 3′ exata; indica possível anelamento, não prova amplificação.',
      'spec.org.label':'Organismo-alvo (busca NCBI)',
      'spec.org.ph':'Digite o nome da espécie (ex.: Babesia bigemina)',
      'btn.listgenomes':'Listar genomas (NCBI)',
      'spec.genome.label':'Genoma / sequência (RefSeq)',
      'spec.genome.none':'— busque um organismo acima —',
      'btn.addgenome':'Adicionar ao fundo',
      'btn.listing':'Listando...',
      'alert.specNoOrg':'Digite ou selecione um organismo (ex.: Babesia bigemina) e tente novamente.',
      'alert.specNoGenome':'Selecione um genoma/sequência na lista (use “Listar genomas (NCBI)” primeiro).',
      'spec.genomeNoneFound':'— nenhuma sequência RefSeq genômica encontrada; tente outro nome ou use um acesso —',
      'spec.genomeOpt': (acc, len, title) => `${acc} · ${len} nt · ${title}`,
      'alert.specGenomeErr':'Falha ao consultar genomas no NCBI (rede). Detalhe: ',
      // varredura de Mg²⁺ → heatmap
      'spec.mg.h3':'Varredura de Mg²⁺ → heatmap',
      'spec.mg.help':'Varre o MgSO₄ total e calcula, para cada primer do Conjunto #1, o Tm pareado (on-target) e o Tm do off-target (melhor janela no fundo). Caminho de Tm independente: Mg livre por quelação do dNTP + correção divalente (Owczarzy 2008) e NN de mismatch interno (Allawi/Peyret). O heatmap colore pela margem de segurança do off-target (T_reação − Tm_off, °C) e sugere o menor [Mg²⁺] que mantém todos os primers funcionais e seguros.',
      'spec.mg.dntp.label':'dNTP total (mM)',
      'spec.mg.mon.label':'Monovalente Na⁺/K⁺ (mM)',
      'spec.mg.trxn.label':'Temperatura de reação (°C)',
      'spec.mg.mgmin.label':'Mg²⁺ mín. (mM)',
      'spec.mg.mgmax.label':'Mg²⁺ máx. (mM)',
      'spec.mg.mgstep.label':'Passo de Mg²⁺ (mM)',
      'spec.mg.run':'Varredura Mg²⁺ → heatmap',
      'spec.mg.primer':'Primer',
      'spec.mg.valid':'Válido',
      'spec.mg.safety':'segurança (T_reação − Tm_off)',
      'spec.mg.status':'status',
      'spec.mg.notFunc':'não-funcional',
      'spec.mg.notSafe':'risco off-target',
      'spec.mg.suggest': (tot, free, minSafety) => `Mg²⁺ sugerido: ${tot} mM total (~${free} mM livre) · margem de segurança mínima do off-target = ${minSafety} °C (todos os primers funcionais e seguros).`,
      'spec.mg.noFeas': (tot, free) => `Nenhuma coluna de Mg²⁺ deixou todos os primers funcionais e seguros. Melhor coluna: ${tot} mM total (~${free} mM livre), maximizando primers funcionais e seguros.`,
      'spec.mg.none':'Sem resultado de varredura. Carregue um fundo e gere primers primeiro.',
      'spec.mg.offenders': names => `Primers funcionais/seguros falhando nesta coluna: ${names}.`,
      'spec.mg.note':'Cor = margem de segurança do off-target (T_reação − Tm_off, °C). Verde = off-target derrete abaixo da T de reação (seguro, maior = mais seguro); vermelho = Tm_off ≥ T de reação (risco de pareamento do off-target); cinza = não-funcional (Tm_on fora da janela LAMP). Anotação da célula = Tm_off absoluto (°C). Caixa azul = Mg²⁺ sugerido.',
      'spec.mg.note2':'Mg desloca o Tm absoluto (rendimento e o Tm do off-target), não a margem de discriminação S; a especificidade é governada pela divergência da sequência.',
    },
    en: {
      'tab.input':'Input', 'tab.params':'Parameters', 'tab.results':'Results',
      'seq.h2':'Target sequence',
      'seq.help':'Paste the sequence (A, T, G, C) or provide an identifier (e.g., GenBank).',
      'id.label':'Identifier (optional)', 'id.ph':'e.g., NC_045512, MN908947, etc.',
      'seq.label':'Sequence', 'seq.ph':'Paste the DNA sequence here...',
      'btn.fasta':'Load FASTA', 'btn.clear':'Clear', 'btn.run':'Generate primers',
      'refs.h3':'References and Publications',
      'refs.lead':'Key references implemented in the LAMPrime algorithm:',
      'refs.tools':'Primer geometry rules: Eiken PrimerExplorer V5 (https://primerexplorer.jp/) and New England Biolabs — LAMP Primer Design (https://lamp.neb.com/).',
      'refs.github':'Open source (MIT): <a href="https://github.com/luizdione/LAMPrime_browser-based" target="_blank" rel="noopener" style="color:var(--accent);">github.com/luizdione/LAMPrime_browser-based</a>',
      'params.basic.h2':'Basic settings',
      'tm.inner.label':'Target Tm (F1c/B1c/LF/LB)', 'tm.outer.label':'Target Tm (F2/B2/F3/B3)',
      'size.f3b3.label':'Primer size (F3/B3)', 'size.lflb.label':'Primer size (LF/LB)',
      'loopseq.label':'Inner loop linker (FIP/BIP) — 5´TTTTTT 3´',
      'fip.label':'Total FIP length (F1c-loop-F2)', 'bip.label':'Total BIP length (B1c-loop-B2)',
      'gc.label':'Expected GC',
      'gc.opt.4060':'40–60%', 'gc.opt.4065':'40–65% (Eiken canonical ceiling)', 'gc.opt.5060':'50–60% (optimal)',
      'gc.opt.3565':'35–65%', 'gc.opt.4555':'45–55%',
      'loops.label':'Loop primers',
      'loops.opt.auto':'Automatic', 'loops.opt.on':'Force inclusion', 'loops.opt.off':'Do not include',
      'lc.label':'Minimum linguistic complexity (LC%)',
      'lc.help':'Linguistic complexity (LC%) measures the “vocabulary richness” of a genetic sequence, from the number of distinct nucleotide combinations (k-mers) observed relative to the theoretical maximum, expressed as a percentage (100% = maximum). Allowed range: 70–90%. Default: 75%.',
      'react.h2':'Reaction conditions',
      'pen.h2':'Penalties and ranking',
      'pen.f1cb1c.label':'F1c/B1c limit', 'pen.f2b2.label':'F2/B2 limit', 'pen.f3b3.label':'F3/B3 limit',
      'maxsets.label':'Maximum number of sets to test',
      'pen.help':'Penalties are used only to rank primer sets. No set is rejected by penalty: the limits are scoring references.',
      'geo.h2':'Target geometry',
      'geo.f2b2.label':'Amplicon F2–B2 (5′F2→5′B2)',
      'geo.opt.120180':'120–180 bp', 'geo.opt.120160':'120–160 bp (PrimerExplorer/NEB)', 'geo.opt.120200':'120–200 bp',
      'geo.f3f2.label':'Gap F3–F2 (3′F3→5′F2)',
      'geo.opt.060':'0–60 bp', 'geo.opt.020':'0–20 bp (Eiken loopamp)', 'geo.opt.040':'0–40 bp',
      'geo.f1b1.label':'F1–B1 distance', 'geo.f2f1.label':'F2–F1 distance', 'geo.b2b1.label':'B2–B1 distance',
      'geo.help':'Geometry per the Eiken PrimerExplorer V5 / NEB standard: <b>F2–B2 amplicon (5′F2→5′B2) = 120–180 bp</b> (amplified region); <b>F3–F2 gap (3′F3→5′F2) = 0–60 bp</b>; F2–F1 and B2–B1 loops = 40–60 bp; F1–B1 center = 0–60 bp.',
      'struct.h2':'Secondary structures (oligo hybridization)',
      'struct.treact':'Reaction temperature', 'struct.hpth':'ΔG threshold — hairpin',
      'struct.dimerth':'ΔG threshold — dimers (homo/hetero)', 'struct.enddg':'End stability (ΔG ≤)',
      'struct.help':'At the reaction temperature, LAMPrime evaluates the structures that hinder LAMP: hairpins (intramolecular) and dimers — oligo hybridization, homo- and hetero-dimers — using the nearest-neighbor model (SantaLucia 1998). Sets with any structure more stable (more negative ΔG) than the threshold are penalized and ranked lower; none are rejected. Priming-end stability follows the Eiken criterion (ΔG ≤ −4 kcal/mol by default).',
      'btn.restore':'Restore defaults', 'btn.apply':'Apply parameters',
      'img.caption':'Schematic of primer arrangement in the LAMP method.',
      'res.h2':'Results',
      'res.empty':'No results yet. Generate primers to see them here.',
      'footer.html':'<strong>LAMPrime</strong> — LAMP primer design in the browser (nearest-neighbor Tm, SantaLucia 1998; PrimerExplorer V5/NEB geometry). Open source (MIT): <a href="https://github.com/luizdione/LAMPrime_browser-based" target="_blank" rel="noopener" style="color:var(--accent);">github.com/luizdione/LAMPrime_browser-based</a>.',
      'footer.lab':'Molecular Genetics Laboratory — IFRJ. De-Melo, Luiz Dione Barbosa — Copyright 2025–2026.',
      // dynamic
      'btn.applying':'Applied!', 'btn.generating':'Generating...',
      'alert.lc':'Enter an LC% value between 70 and 90.',
      'alert.shortSeq':'Sequence too short. Provide at least ~150–250 nt for a reliable LAMP design.',
      'alert.noSets':'No LAMP set satisfied the constraints (Tm/distances). Try a longer sequence (~200–300 nt), or widen the Tm or GC ranges.',
      'alert.designErr':'Design error: ',
      'res.summary':'Design summary', 'res.target':'Target', 'res.sets':'Sets',
      'res.targetTm':'Target Tm', 'res.salts':'Salts',
      'res.method':'Tm by nearest-neighbor (SantaLucia 1998) with salt correction + Mg²⁺ equivalent (von Ahsen 2001). The penalty only ranks; warnings do not reject the set.',
      'res.outer':'Outer primers',
      'res.inner':'Inner primers  (FIP = F1c–loop–F2 · BIP = B1c–loop–B2)',
      'res.loopHdr':'Loop primers (accelerate the reaction)',
      'res.positions':'Positions (1-based)', 'res.distances':'Distances',
      'res.struct': t => `Secondary structure (@${t} °C)`, 'res.worstHp':'worst hairpin ΔG', 'res.worstDi':'worst dimer ΔG',
      'res.allOk':'✓ all primers within Tm/GC/LC ranges', 'res.copy':'copy',
      'res.set':'Set', 'res.penalty':'penalty', 'res.amplicon':'amplicon',
      'warn.hairpin':'hairpin', 'warn.selfdimer':'self-dimer', 'warn.dimer':'dimer',
      'notes.geometry': g => `Geometry (Eiken/NEB standard): F2–B2 amplicon ${g.f2b2[0]}–${g.f2b2[1]} nt · F3–F2 gap ${g.gap[0]}–${g.gap[1]} nt · F2–F1 loop ${g.loopF2F1[0]}–${g.loopF2F1[1]} nt · F1–B1 center ${g.f1b1[0]}–${g.f1b1[1]} nt.`,
      'res.setsOf': n => `(of ${n} valid candidates)`,
      // specificity
      'tab.spec':'Specificity',
      'spec.h2':'Specificity (in-silico screen)',
      'spec.help':'Optional step: screens the primers of the top-ranked set (Set #1) against background sequences, flagging primers that could anneal and prime in related organisms. Uses an exact 3′ seed + extension (k-mer/alignment), not ML.',
      'spec.acc.label':'NCBI accessions (optional)',
      'spec.acc.ph':'e.g., NC_043815.1, AY648786.1 (comma/space)',
      'spec.examples':'Cross-reactivity examples:',
      'spec.examplesHelp':'Fills the congeneric organism\'s NCBI accession; generate primers, fetch the FASTA, then Screen.',
      'btn.fetchncbi':'Fetch FASTA from NCBI',
      'spec.bg.label':'Background sequences (FASTA)',
      'spec.bg.ph':'Paste one or more background FASTA (genomes/genes of related organisms)...',
      'spec.seed.label':'Exact 3′ anchor (seed)',
      'spec.maxmm.label':'Max mismatches (whole primer)',
      'btn.screen':'Screen specificity',
      'spec.privacy':'Privacy: the screen runs in your browser and the primers never leave the page. The network is used only if you click “Fetch FASTA from NCBI” (only the public accession IDs are sent).',
      'spec.res.h2':'Screening result',
      'spec.empty':'No screen yet. Generate primers (Input tab) and provide background sequences.',
      'btn.screening':'Screening...', 'btn.fetching':'Fetching...',
      'alert.specNoDesign':'Generate primers first (Input tab → Generate primers). The screen uses Set #1.',
      'alert.specNoBg':'Provide at least one background sequence (paste FASTA or fetch by NCBI accession).',
      'alert.specNoAcc':'Enter at least one accession number (e.g., NC_043815.1).',
      'alert.specFetchErr':'NCBI fetch failed (CORS or network). Download the FASTA from NCBI manually and paste it. Detail: ',
      'spec.bgHdr': (name, len) => `Background: ${name} (${len} nt)`,
      'spec.setEval': r => `Set screened: #${r}`,
      'spec.summary': (hit, tot) => `${hit} of ${tot} primers with a possible 3′-anchored anneal`,
      'spec.hit': (n, mm, pos, strand) => `⚠ ${n}: ${mm} mismatch(es), pos ${pos} (strand ${strand})`,
      'spec.clean': n => `✓ ${n}: no 3′-anchored anneal`,
      'spec.multi':'Warning: several primers annealing to the same background may indicate cross-amplification risk.',
      'spec.note':'A hit is the best gap-free alignment with an exact 3′ seed; it indicates possible annealing, not proof of amplification.',
      'spec.org.label':'Target organism (NCBI search)',
      'spec.org.ph':'Type the species name (e.g., Babesia bigemina)',
      'btn.listgenomes':'List genomes (NCBI)',
      'spec.genome.label':'Genome / sequence (RefSeq)',
      'spec.genome.none':'— search an organism above —',
      'btn.addgenome':'Add to background',
      'btn.listing':'Listing...',
      'alert.specNoOrg':'Type or select an organism (e.g., Babesia bigemina) and try again.',
      'alert.specNoGenome':'Select a genome/sequence from the list (use "List genomes (NCBI)" first).',
      'spec.genomeNoneFound':'— no RefSeq genomic sequence found; try another name or use an accession —',
      'spec.genomeOpt': (acc, len, title) => `${acc} · ${len} nt · ${title}`,
      'alert.specGenomeErr':'Failed to query genomes at NCBI (network). Detail: ',
      // Mg²⁺ sweep → heatmap
      'spec.mg.h3':'Mg²⁺ sweep → heatmap',
      'spec.mg.help':'Sweeps total MgSO₄ and computes, for each primer of Set #1, the on-target (matched) Tm and the off-target Tm (best window in the background). Independent Tm path: free Mg from dNTP chelation + divalent correction (Owczarzy 2008) and internal-mismatch NN (Allawi/Peyret). The heatmap colors by the off-target safety margin (T_rxn − Tm_off, °C) and suggests the lowest [Mg²⁺] that keeps every primer functional and safe.',
      'spec.mg.dntp.label':'Total dNTP (mM)',
      'spec.mg.mon.label':'Monovalent Na⁺/K⁺ (mM)',
      'spec.mg.trxn.label':'Reaction temperature (°C)',
      'spec.mg.mgmin.label':'Mg²⁺ min (mM)',
      'spec.mg.mgmax.label':'Mg²⁺ max (mM)',
      'spec.mg.mgstep.label':'Mg²⁺ step (mM)',
      'spec.mg.run':'Mg²⁺ sweep → heatmap',
      'spec.mg.primer':'Primer',
      'spec.mg.valid':'Valid',
      'spec.mg.safety':'safety (T_rxn − Tm_off)',
      'spec.mg.status':'status',
      'spec.mg.notFunc':'non-functional',
      'spec.mg.notSafe':'off-target risk',
      'spec.mg.suggest': (tot, free, minSafety) => `Suggested Mg²⁺: ${tot} mM total (~${free} mM free) · minimum off-target safety margin = ${minSafety} °C (all primers functional and safe).`,
      'spec.mg.noFeas': (tot, free) => `No Mg²⁺ column kept all primers functional and safe. Best column: ${tot} mM total (~${free} mM free), maximizing functional and safe primers.`,
      'spec.mg.none':'No sweep result. Load a background and generate primers first.',
      'spec.mg.offenders': names => `Functional/safe primers failing in this column: ${names}.`,
      'spec.mg.note':'Color = off-target safety margin (T_rxn − Tm_off, °C). Green = off-target melts below the reaction temperature (safe, larger = safer); red = Tm_off ≥ reaction temperature (off-target could prime); grey = non-functional (Tm_on outside the LAMP window). Cell annotation = absolute Tm_off (°C). Blue box = suggested Mg²⁺.',
      'spec.mg.note2':'Mg shifts absolute Tm (yield and the off-target\'s Tm), not the match/mismatch margin S; specificity is governed by sequence divergence.',
    }
  };
  let LANG = (function(){ const m=location.search.match(/[?&]lang=(en|pt)/i); if(m) return m[1].toLowerCase(); const hl=(document.documentElement.getAttribute('lang')||'').toLowerCase(); if(hl.indexOf('en')===0) return 'en'; if(hl.indexOf('pt')===0) return 'pt'; try{ const s=localStorage.getItem('lamprime_lang'); if(s==='en'||s==='pt') return s; }catch(e){} return 'pt'; })();
  const L = () => I18N[LANG] || I18N.pt;

  // Tabs
  const tabs = qsa('.lw-tab');
  const panels = qsa('.lw-panel');
  tabs.forEach(tab => {
    tab.addEventListener('click', () => {
      const target = tab.getAttribute('aria-controls');
      tabs.forEach(t => t.classList.toggle('is-active', t === tab));
      tabs.forEach(t => t.setAttribute('aria-selected', t === tab ? 'true' : 'false'));
      panels.forEach(p => p.classList.toggle('is-active', p.id === target));
    });
  });

  // Entrada
  const seqEl = qs('#alvo-seq');
  const idEl = qs('#alvo-id');
  const fileEl = qs('#file-fasta');
  const btnLimpar = qs('#btn-limpar');
  const btnExecutar = qs('#btn-executar');

  btnLimpar.addEventListener('click', () => {
    seqEl.value = '';
    idEl.value = '';
  });

  fileEl.addEventListener('change', async (e) => {
    const file = e.target.files && e.target.files[0];
    if (!file) return;
    const text = await file.text();
    // FASTA simples: remove headers '>' e concatena linhas
    const seq = text.split(/\r?\n/).filter(l => !l.startsWith('>')).join('').trim();
    seqEl.value = seq;
  });

  // Parâmetros
  const defaults = {
    tmF1cB1cLfLbMin: 64,
    tmF1cB1cLfLbMax: 66,
    tmF2B2F3B3Min: 59,
    tmF2B2F3B3Max: 61,
    tamF3B3Min: 18,
    tamF3B3Max: 23,
    tamLfLbMin: 18,
    tamLfLbMax: 23,
    loopSeq: 'TTTTTT',
    tamFipMin: 42,
    tamFipMax: 52,
    tamBipMin: 42,
    tamBipMax: 52,
    gcRange: '40-60',
    loops: 'auto',
    lcMin: 75,
    na: 50,
    mg: 8,
    dntp: 1.4,
    distF3F2: '0-60',
    distF2B2: '120-180',
    distF1B1: '0-60',
    distF2F1: '40-60',
    distB2B1: '40-60',
    penLimF1cB1c: 15,
    penLimF2B2: 20,
    penLimF3B3: 15,
    maxSets: 300,
    tReact: 63,
    hpTh: -3.0,
    dimerTh: -3.0,
    endDg: -4.0
  };

  const els = {
    tmF1cB1cLfLbMin: qs('#tm-f1c-b1c-lf-lb-min'),
    tmF1cB1cLfLbMax: qs('#tm-f1c-b1c-lf-lb-max'),
    tmF2B2F3B3Min: qs('#tm-f2-b2-f3-b3-min'),
    tmF2B2F3B3Max: qs('#tm-f2-b2-f3-b3-max'),
    tamF3B3Min: qs('#tam-f3-b3-min'),
    tamF3B3Max: qs('#tam-f3-b3-max'),
    tamLfLbMin: qs('#tam-lf-lb-min'),
    tamLfLbMax: qs('#tam-lf-lb-max'),
    loopSeq: qs('#loop-seq'),
    tamFipMin: qs('#tam-fip-min'),
    tamFipMax: qs('#tam-fip-max'),
    tamBipMin: qs('#tam-bip-min'),
    tamBipMax: qs('#tam-bip-max'),
    gcRange: qs('#gc-range'),
    loops: qs('#loops'),
    lcMin: qs('#lc-min'),
    na: qs('#na'),
    mg: qs('#mg'),
    dntp: qs('#dnTP'),
    distF3F2: qs('#dist-f3-f2'),
    distF2B2: qs('#dist-f2-b2'),
    distF1B1: qs('#dist-f1-b1'),
    distF2F1: qs('#dist-f2-f1'),
    distB2B1: qs('#dist-b2-b1'),
    penLimF1cB1c: qs('#pen-lim-f1c-b1c'),
    penLimF2B2: qs('#pen-lim-f2-b2'),
    penLimF3B3: qs('#pen-lim-f3-b3'),
    maxSets: qs('#max-sets'),
    tReact: qs('#t-react'),
    hpTh: qs('#hp-th'),
    dimerTh: qs('#dimer-th'),
    endDg: qs('#end-dg'),
  };

  const btnPadrao = qs('#btn-param-padrao');
  const btnAplicar = qs('#btn-aplicar-param');

  function restoreDefaults() {
    Object.entries(defaults).forEach(([k, v]) => {
      if (els[k]) {
        els[k].value = v;
        const ev = new Event('change', { bubbles: true });
        els[k].dispatchEvent(ev);
      }
    });
  }

  btnPadrao.addEventListener('click', restoreDefaults);
  // Normaliza sequência do loop para letras maiúsculas
  if (els.loopSeq) {
    els.loopSeq.addEventListener('input', () => {
      els.loopSeq.value = (els.loopSeq.value || '').toUpperCase().replace(/[^ATGC]/g, '');
    });
  }
  btnAplicar.addEventListener('click', () => {
    // Validação simples do LC%
    const v = parseFloat(els.lcMin.value);
    if (isNaN(v) || v < 70 || v > 90) {
      alert(L()['alert.lc']);
      els.lcMin.focus();
      return;
    }
    // Feedback simples
    btnAplicar.textContent = L()['btn.applying'];
    setTimeout(() => btnAplicar.textContent = L()['btn.apply'], 1000);
  });

  restoreDefaults();

  // Resultados
  const resEmpty = qs('#res-vazio');
  const resList = qs('#res-list');
  let lastData = null; // último resultado renderizado (para re-render ao trocar idioma)

  function toTabResultados() {
    const tabRes = qs('#tabbtn-res');
    tabRes.click();
  }

  function sanitizeDNA(seq) {
    return (seq || '').toUpperCase().replace(/U/g, 'T').replace(/[^ATGC]/g, '');
  }

  // ===================================================================
  // Bioinformática — Tm nearest-neighbor (SantaLucia 1998) + sais
  // Refs: Notomi 2000; Nagamine 2002; Tomita 2008; SantaLucia 1998;
  //       von Ahsen 2001 (eq. Mg→Na+); PrimerExplorer V5 (Eiken).
  // ===================================================================
  const COMP = { A:'T', T:'A', G:'C', C:'G', N:'N' };
  const comp = b => COMP[b] || 'N';
  const revComp = s => { let o=''; for (let i=s.length-1;i>=0;i--) o+=comp(s[i]); return o; };
  const gcPct = s => { const n=s.length||1; return 100*((s.match(/[GC]/g)||[]).length)/n; };

  // ΔH (kcal/mol) e ΔS (cal/mol·K) — pares vizinhos 5'XY3' (SantaLucia 1998 unificado)
  const NN_DH = {AA:-7.9,AT:-7.2,TA:-7.2,CA:-8.5,GT:-8.4,CT:-7.8,GA:-8.2,CG:-10.6,GC:-9.8,GG:-8.0,AC:-8.4,AG:-7.8,TC:-8.2,TG:-8.5,TT:-7.9,CC:-8.0};
  const NN_DS = {AA:-22.2,AT:-20.4,TA:-21.3,CA:-22.7,GT:-22.4,CT:-21.0,GA:-22.2,CG:-27.2,GC:-24.4,GG:-19.9,AC:-22.4,AG:-21.0,TC:-22.2,TG:-22.7,TT:-22.2,CC:-19.9};
  const R_GAS = 1.987;
  const TM_OLIGO_NM = 50; // conc. representativa para estimativa de Tm (ajustável)

  // Mg²⁺ livre via equilíbrio de quelação 1:1 com dNTP (Ka=3e4 M⁻¹). Conc. em MOLAR.
  // Helper COMPARTILHADO (#6): usado tanto pelo motor de desenho (tmNN, via von Ahsen 2001)
  // quanto pela varredura de especificidade (tmDuplexMg, via Owczarzy 2008) — unifica o Mg livre.
  function freeMgM(totMgM, dntpM, Ka=3e4){ const B=Ka*dntpM - Ka*totMgM + 1; return (-(B)+Math.sqrt(B*B+4*Ka*totMgM))/(2*Ka); }

  function tmNN(seq, naMM, mgMM, dntpMM, oligoNM) {
    seq = (seq||'').toUpperCase().replace(/[^ATGC]/g,'');
    const N = seq.length; if (N < 2) return NaN;
    let dH=0, dS=0;
    for (let i=0;i<N-1;i++){ const d=seq.substr(i,2); if(NN_DH[d]===undefined) return NaN; dH+=NN_DH[d]; dS+=NN_DS[d]; }
    const init = b => (b==='G'||b==='C') ? [0.1,-2.8] : [2.3,4.1];
    const [h5,s5]=init(seq[0]); const [h3,s3]=init(seq[N-1]); dH+=h5+h3; dS+=s5+s3;
    // Na+ equivalente incluindo Mg2+ livre (von Ahsen 2001): Na_eq = Na + 120*sqrt(Mg_livre).
    // #6: Mg2+ livre pela quelação de equilíbrio (freeMgM, Ka=3e4), unificado com a varredura,
    // no lugar da subtração crua max(0, Mg-dNTP) — que zerava o Mg livre na estequiometria.
    const mgFree = freeMgM((mgMM||0)/1000, (dntpMM||0)/1000) * 1000; // M -> mM
    const naEq = Math.max(1e-3, (naMM||0) + 120*Math.sqrt(Math.max(0,mgFree))) / 1000; // mol/L
    const dS_salt = dS + 0.368*(N-1)*Math.log(naEq); // correção de sal SantaLucia 1998
    const C = (oligoNM||TM_OLIGO_NM)*1e-9;
    return (dH*1000) / (dS_salt + R_GAS*Math.log(C/4)) - 273.15;
  }

  // Complexidade linguística LC% — média geométrica da diversidade de k-mers (k=1..3)
  function linguisticComplexity(s) {
    s=(s||'').toUpperCase(); const N=s.length; if (N<3) return 100;
    let prod=1, terms=0;
    for (let k=1;k<=3;k++){
      const max=Math.min(Math.pow(4,k), N-k+1); if (max<=0) continue;
      const set=new Set(); for (let i=0;i+k<=N;i++) set.add(s.substr(i,k));
      prod *= set.size/max; terms++;
    }
    return 100*Math.pow(prod, terms?1/terms:1);
  }

  // ΔG (37°C) dos ~6 nt da extremidade — critério Eiken: extremidades de iniciação
  // devem ser estáveis (ΔG ≤ -4 kcal/mol). end: '3' (3') ou '5' (5').
  function endDeltaG(seq, end) {
    seq=(seq||'').toUpperCase().replace(/[^ATGC]/g,'');
    const w = (end==='5') ? seq.slice(0,6) : seq.slice(-6);
    if (w.length<2) return 0;
    let dH=0,dS=0;
    for (let i=0;i<w.length-1;i++){ const d=w.substr(i,2); if(NN_DH[d]===undefined) continue; dH+=NN_DH[d]; dS+=NN_DS[d]; }
    return dH - 310.15*(dS/1000); // kcal/mol
  }

  // ===== Estruturas secundárias: hairpin, homodímero, heterodímero =====
  // ΔG (kcal/mol) na TEMPERATURA DE REAÇÃO (~63 °C): só estruturas estáveis nessa
  // faixa (60–65 °C) atrapalham a LAMP. Modelo NN SantaLucia 1998.
  let REACT_C = 63;                         // temperatura de reação (°C), ajustável pelo usuário
  const TREACT_K = () => (REACT_C + 273.15);
  function segDG(seg) { // ΔG de um duplexo perfeitamente pareado de `seg`
    seg = (seg||'').toUpperCase(); const N = seg.length; if (N < 2) return 0;
    let dH=0, dS=0;
    for (let i=0;i<N-1;i++){ const d=seg.substr(i,2); if(NN_DH[d]===undefined) continue; dH+=NN_DH[d]; dS+=NN_DS[d]; }
    const ini = b => (b==='G'||b==='C') ? [0.1,-2.8] : [2.3,4.1];
    const [h5,s5]=ini(seg[0]); const [h3,s3]=ini(seg[N-1]); dH+=h5+h3; dS+=s5+s3;
    return dH - TREACT_K()*dS/1000;
  }
  // duplexo antiparalelo mais estável entre a e b (b=a → homodímero)
  function dimerDG(a, b) {
    a=(a||'').toUpperCase(); const B=(b||'').toUpperCase().split('').reverse().join('');
    let best=0;
    for (let off=-(a.length-1); off<B.length; off++){
      let i=Math.max(0,-off);
      while (i<a.length && i+off<B.length){
        if (i+off>=0 && a[i]===comp(B[i+off])){
          let j=i, seg='';
          while (j<a.length && j+off>=0 && j+off<B.length && a[j]===comp(B[j+off])){ seg+=a[j]; j++; }
          if (seg.length>=3){ const dg=segDG(seg); if (dg<best) best=dg; }
          i=j;
        } else i++;
      }
    }
    return best;
  }
  // hairpin mais estável (haste ≥3 pb, alça ≥3 nt)
  function hairpinDG(a) {
    a=(a||'').toUpperCase(); let best=0;
    for (let i=0;i<a.length;i++){
      for (let j=a.length-1;j>i+3;j--){
        if (a[i]!==comp(a[j])) continue;
        let k=0, seg='';
        while (i+k < j-k && (j-k)-(i+k)-1>=3 && a[i+k]===comp(a[j-k])){ seg+=a[i+k]; k++; }
        if (seg.length>=3){ const dg=segDG(seg); if (dg<best) best=dg; }
      }
    }
    return best;
  }
  // varre estruturas dos primers do conjunto: hairpins, self e cross-dímeros.
  // Avisos retornados como tokens estruturados (traduzidos no render).
  // limiares de ΔG (kcal/mol) na temperatura de reação: mais negativo = problemático.
  // hpTh (hairpin intramolecular) e diTh (dímeros homo/hetero) vêm dos parâmetros.
  function structureScan(primers, hpTh, diTh) {
    let pen=0; const warns=[]; let worstHp=0, worstDi=0;
    for (const p of primers){
      const hp=hairpinDG(p.seq); if (hp<worstHp) worstHp=hp;
      if (hp<hpTh){ pen += (hpTh-hp); warns.push({t:'hairpin', n:p.name, v:hp.toFixed(1)}); }
      const sd=dimerDG(p.seq, p.seq); if (sd<worstDi) worstDi=sd;
      if (sd<diTh){ pen += (diTh-sd); warns.push({t:'selfdimer', n:p.name, v:sd.toFixed(1)}); }
    }
    for (let a=0;a<primers.length;a++) for (let b=a+1;b<primers.length;b++){
      const cd=dimerDG(primers[a].seq, primers[b].seq); if (cd<worstDi) worstDi=cd;
      if (cd<diTh){ pen += (diTh-cd)*0.8; warns.push({t:'dimer', n:`${primers[a].name}/${primers[b].name}`, v:cd.toFixed(1)}); }
    }
    return { penalty: pen*3, warns, worstHairpin: worstHp, worstDimer: worstDi };
  }

  // ===================================================================
  // Mg²⁺-vs-especificidade — caminho de Tm SEPARADO (NÃO usa tmNN). Calcula
  // a margem de Tm on-target − off-target de cada primer ao longo de uma
  // varredura de MgSO₄ total. Usa Mg livre por quelação 1:1 do dNTP +
  // correção divalente de Owczarzy 2008 e NN de mismatch interno
  // (Allawi&SantaLucia 1997/1998; Peyret 1999). ADITIVO: não toca no
  // caminho de desenho/concordância nem nos termos de iniciação publicados.
  // ===================================================================
  // NN de mismatch interno único (Allawi&SantaLucia 1997 G·T; 1998 G·A,C·T[NAR 26:2694],A·C; Peyret 1999 A·A/C·C/G·G/T·T). [dH kcal/mol, dS cal/mol·K], 1 M NaCl. Chave 'XY/WZ' = 5'-XY-3' sobre 3'-WZ-5'.
  const MM_IMM = {
'AG/TT':[1.0,0.9],'AT/TG':[-2.5,-8.3],'CG/GT':[-4.1,-11.7],'CT/GG':[-2.8,-8.0],'GG/CT':[3.3,10.4],'GG/TT':[5.8,16.3],'GT/CG':[-4.4,-12.3],'GT/TG':[4.1,9.5],'TG/AT':[-0.1,-1.7],'TG/GT':[-1.4,-6.2],'TT/AG':[-1.3,-5.3],
'AA/TG':[-0.6,-2.3],'AG/TA':[-0.7,-2.3],'CA/GG':[-0.7,-2.3],'CG/GA':[-4.0,-13.2],'GA/CG':[-0.6,-1.0],'GG/CA':[0.5,3.2],'TA/AG':[0.7,0.7],'TG/AA':[3.0,7.4],
'AC/TT':[0.7,0.2],'AT/TC':[-1.2,-6.2],'CC/GT':[-0.8,-4.5],'CT/GC':[-1.5,-6.1],'GC/CT':[2.3,5.4],'GT/CC':[5.2,13.5],'TC/AT':[1.2,0.7],'TT/AC':[1.0,0.7],
'AA/TC':[2.3,4.6],'AC/TA':[5.3,14.6],'CA/GC':[1.9,3.7],'CC/GA':[0.6,-0.6],'GA/CC':[5.2,14.2],'GC/CA':[-0.7,-3.8],'TA/AC':[3.4,8.0],'TC/AA':[7.6,20.2],
'AA/TA':[1.2,1.7],'CA/GA':[-0.9,-4.2],'GA/CA':[-2.9,-9.8],'TA/AA':[4.7,12.9],
'AC/TC':[0.0,-4.4],'CC/GC':[-1.5,-7.2],'GC/CC':[3.6,8.9],'TC/AC':[6.1,16.4],
'AG/TG':[-3.1,-9.5],'CG/GG':[-4.9,-15.3],'GG/CG':[-6.0,-15.8],'TG/AG':[1.6,3.6],
'AT/TT':[-2.7,-10.8],'CT/GT':[-5.0,-15.8],'GT/CT':[-2.2,-8.4],'TT/AT':[0.2,-1.5]
  };

  // freeMgM (quelação 1:1, Ka=3e4) foi promovida a helper compartilhado junto a tmNN (#6);
  // aqui é usada pela correção divalente de Owczarzy 2008.
  // Correção divalente Owczarzy 2008. tmK1M = Tm(K) a 1 M Na+ (sem termo de sal). monM=[Na+K+Tris/2]. Retorna Tm em °C.
  function owczarzy2008(tmK1M, freeMg, monM, fGC, Nbp){
    const ln=Math.log; let a=3.92,b=-0.911,c=6.26,d=1.42,e=-48.2,f=52.5,g=8.31;
    const R = Math.sqrt(Math.max(freeMg,1e-9))/Math.max(monM,1e-9);
    if (R < 0.22){ /* monovalente dominante: sem termo divalente (LAMP raramente aqui) */ return tmK1M-273.15 + 0; }
    if (R < 6.0){ a=3.92*(0.843-0.352*Math.sqrt(monM)*ln(monM)); d=1.42*(1.279-4.03e-3*ln(monM)-8.03e-3*ln(monM)**2); g=8.31*(0.486-0.258*ln(monM)+5.25e-3*ln(monM)**3); }
    const lm=ln(freeMg);
    const corr=(a + b*lm + fGC*(c + d*lm) + (1/(2*(Nbp-1)))*(e + f*lm + g*lm*lm))*1e-5;
    const tmKc = 1/(1/tmK1M + corr);
    return tmKc - 273.15;
  }
  // Tm de um duplexo `top` (5'->3') pareado com `bottom` (3'->5', i.e. o alvo alinhado lido 3'->5'),
  // posição i pareia top[i]:bottom[i]; permite mismatches internos. opt={mg,dntp,mon,oligoNM}, conc. em mM (exceto oligoNM em nM).
  // Retorna {tm:°C, flag:bool} — flag=true se algum stack não-WC ficou sem parâmetro (penalizado).
  function tmDuplexMg(top, bottom, opt){
    top=(top||'').toUpperCase(); bottom=(bottom||'').toUpperCase();
    const N=Math.min(top.length, bottom.length); if (N<2) return {tm:NaN, flag:true};
    const isWC = (x,y) => (x==='A'&&y==='T')||(x==='T'&&y==='A')||(x==='G'&&y==='C')||(x==='C'&&y==='G');
    let dH=0, dS=0, flag=false, gcMatch=0;
    for (let i=0;i<N;i++){ if (isWC(top[i], bottom[i])) gcMatch += (top[i]==='G'||top[i]==='C')?1:0; }
    for (let i=0;i<N-1;i++){
      const a=top[i], b=top[i+1], w=bottom[i], z=bottom[i+1];
      if (isWC(a,w) && isWC(b,z)){
        // stack Watson-Crick: usa a tabela NN existente (5'XY3' do top)
        const d=a+b; if (NN_DH[d]===undefined){ flag=true; continue; }
        dH+=NN_DH[d]; dS+=NN_DS[d];
      } else {
        // stack com pelo menos um mismatch interno: MM_IMM
        const key=a+b+'/'+w+z;
        let p=MM_IMM[key];
        if (!p){ // tenta a equivalência por inversão de orientação (lê o stack pelo outro sentido)
          const rkey=z+w+'/'+b+a; p=MM_IMM[rkey];
        }
        if (p){ dH+=p[0]; dS+=p[1]; }
        else { // par desconhecido (ex.: mismatch terminal ou duplo mismatch): fortemente desestabilizante
          dH+=0; dS+=12; flag=true; // ΔS positivo penaliza (reduz |dS_total|) → baixa a Tm
        }
      }
    }
    // termos de iniciação (mesmos valores publicados) pelos pares terminais; só pares WC contam
    const ini = (x,y) => isWC(x,y) ? ((x==='G'||x==='C')?[0.1,-2.8]:[2.3,4.1]) : [0,0];
    const [h5,s5]=ini(top[0],bottom[0]); const [h3,s3]=ini(top[N-1],bottom[N-1]); dH+=h5+h3; dS+=s5+s3;
    const C=((opt&&opt.oligoNM)||TM_OLIGO_NM)*1e-9;
    const denom = dS + R_GAS*Math.log(C/4);
    if (!(denom<0)) return {tm:NaN, flag:true};
    const tmK1M = (dH*1000)/denom;            // Tm(K) a 1 M Na+ (sem termo de sal)
    const totMgM=((opt&&opt.mg!=null?opt.mg:8))/1000;  // mM -> M
    const dntpM =((opt&&opt.dntp!=null?opt.dntp:1.4))/1000;
    const monM  =Math.max(1e-3, ((opt&&opt.mon!=null?opt.mon:50))/1000);
    const fGC = gcMatch / N;
    const tmC = owczarzy2008(tmK1M, freeMgM(totMgM, dntpM), monM, fGC, N);
    return {tm:tmC, flag};
  }
  // complemento Watson-Crick alinhado 3'->5' (para a posição i parear top[i]) — alvo perfeito on-target
  function wcBottom(top){ top=(top||'').toUpperCase(); let o=''; for (let i=0;i<top.length;i++) o+=comp(top[i]); return o; }

  // Melhor janela off-target sem gaps: desliza `primer` sobre cada fundo (ambas as fitas),
  // acha a janela com MENOS mismatches; retorna a base do alvo pareada (3'->5') p/ tmDuplexMg.
  // Força bruta O(primerLen*bgLen) — ok p/ amplicon/gene/plasmídeo. (genoma inteiro pediria índice k-mer; trabalho futuro.)
  function bestOffWindow(primer, bgSeqs){
    primer=(primer||'').toUpperCase(); const Lp=primer.length; if (Lp<2) return null;
    let best=null;
    const scan=(seq, strand)=>{
      for (let s=0; s+Lp<=seq.length; s++){
        let mm=0; for (let i=0;i<Lp;i++){ if (seq[s+i]!==primer[i]) mm++; }
        if (!best || mm<best.nMismatch){
          // bottom alinhado 3'->5' tal que top[i]:bottom[i] são as bases efetivamente justapostas:
          // o alvo (seq) fornece a base complementar; bottom[i] = complemento da base do alvo na fita oposta
          let bottom=''; for (let i=0;i<Lp;i++) bottom+=comp(seq[s+i]);
          best={ alignedTargetBottom:bottom, nMismatch:mm, pos:s, strand };
          if (mm===0) return; // ótimo perfeito, pode parar esta fita
        }
      }
    };
    (bgSeqs||[]).forEach(bg=>{
      const fwd=bg.seq||bg; scan(fwd, '+');
      scan(revComp(fwd), '−');
    });
    return best;
  }

  const inRange = (x, lo, hi) => x>=lo && x<=hi;
  const parseRange = (str, dlo, dhi) => { const m=String(str||'').match(/(-?\d+(?:\.\d+)?)\s*-\s*(-?\d+(?:\.\d+)?)/); return m?[parseFloat(m[1]),parseFloat(m[2])]:[dlo,dhi]; };

  // melhor primer começando em `start` (5' na fita sense), comprimento [lenMin,lenMax]
  function bestByStart(S, start, lenMin, lenMax, tmLo, tmHi, P) {
    const mid=(tmLo+tmHi)/2; let best=null;
    for (let L=lenMin; L<=lenMax; L++){
      if (start+L>S.length) break;
      const seg=S.substr(start,L); const tm=tmNN(seg,P.na,P.mg,P.dntp,TM_OLIGO_NM); if (isNaN(tm)) continue;
      const dev=Math.abs(tm-mid);
      if (!best || dev<best.dev) best={start,end:start+L-1,len:L,seg,tm,gc:gcPct(seg),lc:linguisticComplexity(seg),dev};
    }
    return best;
  }
  // melhor primer terminando em `end`
  function bestByEnd(S, end, lenMin, lenMax, tmLo, tmHi, P) {
    const mid=(tmLo+tmHi)/2; let best=null;
    for (let L=lenMin; L<=lenMax; L++){
      const start=end-L+1; if (start<0) break;
      const seg=S.substr(start,L); const tm=tmNN(seg,P.na,P.mg,P.dntp,TM_OLIGO_NM); if (isNaN(tm)) continue;
      const dev=Math.abs(tm-mid);
      if (!best || dev<best.dev) best={start,end,len:L,seg,tm,gc:gcPct(seg),lc:linguisticComplexity(seg),dev};
    }
    return best;
  }

  function readParams() {
    const num=(el,d)=>{const v=parseFloat(el&&el.value); return isNaN(v)?d:v;};
    const [gcLo,gcHi]=parseRange(els.gcRange&&els.gcRange.value,40,60);
    return {
      innerLo:num(els.tmF1cB1cLfLbMin,64), innerHi:num(els.tmF1cB1cLfLbMax,66),
      outerLo:num(els.tmF2B2F3B3Min,59), outerHi:num(els.tmF2B2F3B3Max,61),
      f3Min:num(els.tamF3B3Min,18), f3Max:num(els.tamF3B3Max,23),
      lfMin:num(els.tamLfLbMin,18), lfMax:num(els.tamLfLbMax,23),
      fipMin:num(els.tamFipMin,42), fipMax:num(els.tamFipMax,52),
      bipMin:num(els.tamBipMin,42), bipMax:num(els.tamBipMax,52),
      gcLo, gcHi,
      loops:(els.loops&&els.loops.value)||'auto',
      lcMin:num(els.lcMin,75),
      na:num(els.na,50), mg:num(els.mg,8), dntp:num(els.dntp,1.4),
      d_F2F1:parseRange(els.distF2F1&&els.distF2F1.value,40,60),
      d_F1B1:parseRange(els.distF1B1&&els.distF1B1.value,0,60),
      d_F2B2:parseRange(els.distF2B2&&els.distF2B2.value,120,180),  // amplicon 5'F2→5'B2 (Eiken/NEB: 120–180)
      d_B2B1:parseRange(els.distB2B1&&els.distB2B1.value,40,60),
      gapF3F2:parseRange(els.distF3F2&&els.distF3F2.value,0,60),    // gap real 3'F3→5'F2 (Eiken/NEB: 0–60)
      loopLinker:((els.loopSeq&&els.loopSeq.value)||'').toUpperCase().replace(/[^ATGC]/g,''),
      limF1cB1c:num(els.penLimF1cB1c,15), limF2B2:num(els.penLimF2B2,20), limF3B3:num(els.penLimF3B3,15),
      maxSets:Math.max(1,Math.min(10000,parseInt((els.maxSets&&els.maxSets.value)||300,10))),
      tReact:num(els.tReact,63), hpTh:num(els.hpTh,-3.0), dimerTh:num(els.dimerTh,-3.0), endDg:num(els.endDg,-4.0)
    };
  }

  function designLAMP(S, P) {
    const N=S.length;
    if (N<120) return {errorKey:'alert.shortSeq'};
    REACT_C = (P.tReact!=null) ? P.tReact : 63; // ΔG das estruturas 2ª na temperatura de reação
    const innerLenMin=18, innerLenMax=26;
    const geom={ f2b2:P.d_F2B2, gap:P.gapF3F2, loopF2F1:P.d_F2F1, f1b1:P.d_F1B1 };

    const blockDev = b => (b.F1?b.F1.dev:b.B1.dev)+(b.F2?b.F2.dev:b.B2.dev)+(b.F3?b.F3.dev:b.B3.dev);

    // 1) F-blocks: F3 < F2 < F1 (sentido 5'->3')
    let Fblocks=[];
    for (let f2s=0; f2s<N; f2s++) {
      const F2=bestByStart(S,f2s,P.f3Min,P.f3Max,P.outerLo,P.outerHi,P);
      if (!F2 || !inRange(F2.tm,P.outerLo,P.outerHi)) continue;
      let F1=null;
      for (let d=P.d_F2F1[0]; d<=P.d_F2F1[1]; d++){
        const c=bestByStart(S,f2s+d,innerLenMin,innerLenMax,P.innerLo,P.innerHi,P);
        if (c && inRange(c.tm,P.innerLo,P.innerHi) && (!F1||c.dev<F1.dev)) F1=c;
      }
      if (!F1) continue;
      let F3=null;
      for (let g=P.gapF3F2[0]; g<=P.gapF3F2[1]; g++){
        const f3e=f2s-1-g; if (f3e<P.f3Min-1) break;
        const c=bestByEnd(S,f3e,P.f3Min,P.f3Max,P.outerLo,P.outerHi,P);
        if (c && inRange(c.tm,P.outerLo,P.outerHi) && (!F3||c.dev<F3.dev)) F3=c;
      }
      if (!F3) continue;
      Fblocks.push({F3,F2,F1});
    }
    Fblocks.sort((a,b)=>blockDev(a)-blockDev(b)); Fblocks=Fblocks.slice(0,250);

    // 2) B-blocks na fita sense: B1 < B2 < B3; B1 inner, B2/B3 outer
    let Bblocks=[];
    for (let b1s=0; b1s<N; b1s++) {
      const B1=bestByStart(S,b1s,innerLenMin,innerLenMax,P.innerLo,P.innerHi,P);
      if (!B1 || !inRange(B1.tm,P.innerLo,P.innerHi)) continue;
      let B2=null;
      for (let d=P.d_B2B1[0]; d<=P.d_B2B1[1]; d++){
        const c=bestByStart(S,b1s+d,P.f3Min,P.f3Max,P.outerLo,P.outerHi,P);
        if (c && inRange(c.tm,P.outerLo,P.outerHi) && (!B2||c.dev<B2.dev)) B2=c;
      }
      if (!B2) continue;
      let B3=null;
      for (let g=P.gapF3F2[0]; g<=P.gapF3F2[1]; g++){
        const c=bestByStart(S,B2.end+1+g,P.f3Min,P.f3Max,P.outerLo,P.outerHi,P);
        if (c && inRange(c.tm,P.outerLo,P.outerHi) && (!B3||c.dev<B3.dev)) B3=c;
      }
      if (!B3) continue;
      Bblocks.push({B1,B2,B3});
    }
    Bblocks.sort((a,b)=>blockDev(a)-blockDev(b)); Bblocks=Bblocks.slice(0,250);

    // 3) combinar respeitando F1-B1 (centro) e F2-B2 (amplicon)
    const sets=[];
    const midI=(P.innerLo+P.innerHi)/2, midO=(P.outerLo+P.outerHi)/2;
    for (const F of Fblocks){
      for (const B of Bblocks){
        if (!(F.F1.end < B.B1.start)) continue;            // F1 à esquerda de B1
        const dF1B1 = B.B1.start - F.F1.end - 1;
        if (!inRange(dF1B1, P.d_F1B1[0], P.d_F1B1[1])) continue;
        const dF2B2 = B.B2.end - F.F2.start + 1;
        if (!inRange(dF2B2, P.d_F2B2[0], P.d_F2B2[1])) continue;
        if (!(F.F3.start<F.F2.start && F.F2.end<F.F1.start && B.B1.end<B.B2.start && B.B2.end<B.B3.start)) continue;

        // loops
        const lfRegion = S.substring(F.F2.end+1, F.F1.start);
        const lbRegion = S.substring(B.B1.end+1, B.B2.start);
        let LF=null, LB=null;
        if (P.loops!=='off'){
          if (lfRegion.length>=P.lfMin){ const c=bestByStart(lfRegion,0,P.lfMin,Math.min(P.lfMax,lfRegion.length),P.innerLo,P.innerHi,P); if(c) LF={seq:revComp(c.seg),tm:c.tm,gc:c.gc,len:c.len}; }
          if (lbRegion.length>=P.lfMin){ const c=bestByStart(lbRegion,0,P.lfMin,Math.min(P.lfMax,lbRegion.length),P.innerLo,P.innerHi,P); if(c) LB={seq:c.seg,tm:c.tm,gc:c.gc,len:c.len}; }
        }

        const F1c=revComp(F.F1.seg), F2=F.F2.seg, B1c=B.B1.seg, B2=revComp(B.B2.seg);
        const FIP=F1c+P.loopLinker+F2, BIP=B1c+P.loopLinker+B2;
        const F3=F.F3.seg, B3=revComp(B.B3.seg);

        const pInner=(Math.abs(F.F1.tm-midI)+Math.abs(B.B1.tm-midI))/Math.max(1e-6,P.limF1cB1c);
        const pF2B2=(Math.abs(F.F2.tm-midO)+Math.abs(B.B2.tm-midO))/Math.max(1e-6,P.limF2B2);
        const pF3B3=(Math.abs(F.F3.tm-midO)+Math.abs(B.B3.tm-midO))/Math.max(1e-6,P.limF3B3);
        const segs=[['F3',F.F3],['F2',F.F2],['F1c',F.F1],['B1c',B.B1],['B2',B.B2],['B3',B.B3]];
        let gcPen=0,lcPen=0,warns=[];
        for (const [name,info] of segs){
          if (!inRange(info.gc,P.gcLo,P.gcHi)){ gcPen+=1; warns.push({t:'gc', n:name, v:info.gc.toFixed(0)}); }
          if (info.lc<P.lcMin){ lcPen+=1; warns.push({t:'lc', n:name, v:info.lc.toFixed(0)}); }
        }
        const fipOk=inRange(FIP.length,P.fipMin,P.fipMax), bipOk=inRange(BIP.length,P.bipMin,P.bipMax);
        if (!fipOk) warns.push({t:'len', n:'FIP', v:FIP.length}); if (!bipOk) warns.push({t:'len', n:'BIP', v:BIP.length});
        // estabilidade das extremidades de iniciação (Eiken: ΔG ≤ -4 kcal/mol)
        let dgPen=0;
        [['F3',F3,'3'],['F2',F2,'3'],['B3',B3,'3'],['B2',B2,'3']].forEach(([n,s,e])=>{ const g=endDeltaG(s,e); if(g>P.endDg){dgPen++; warns.push({t:'dg', n, end:e, v:g.toFixed(1)});} });
        [['F1c',F1c,'5'],['B1c',B1c,'5']].forEach(([n,s,e])=>{ const g=endDeltaG(s,e); if(g>P.endDg){dgPen++; warns.push({t:'dg', n, end:e, v:g.toFixed(1)});} });
        const penalty=(pInner+pF2B2+pF3B3)*100 + gcPen*5 + lcPen*5 + dgPen*4 + (fipOk?0:8) + (bipOk?0:8);

        sets.push({
          penalty,
          ampliconF3B3: B.B3.end - F.F3.start + 1, ampliconF2B2: dF2B2, dF1B1,
          coords:{F3:[F.F3.start+1,F.F3.end+1],F2:[F.F2.start+1,F.F2.end+1],F1:[F.F1.start+1,F.F1.end+1],B1:[B.B1.start+1,B.B1.end+1],B2:[B.B2.start+1,B.B2.end+1],B3:[B.B3.start+1,B.B3.end+1]},
          F3:{seq:F3,tm:F.F3.tm,gc:F.F3.gc}, B3:{seq:B3,tm:B.B3.tm,gc:B.B3.gc},
          F2:{seq:F2,tm:F.F2.tm,gc:F.F2.gc}, B2:{seq:B2,tm:B.B2.tm,gc:B.B2.gc},
          F1c:{seq:F1c,tm:F.F1.tm,gc:F.F1.gc}, B1c:{seq:B1c,tm:B.B1.tm,gc:B.B1.gc},
          FIP:{seq:FIP,len:FIP.length}, BIP:{seq:BIP,len:BIP.length},
          LF, LB, warns
        });
      }
    }
    sets.sort((a,b)=>a.penalty-b.penalty);
    const seen=new Set(), uniq=[];
    for (const s of sets){ const k=s.coords.F2[0]+'|'+s.coords.B2[1]; if(seen.has(k))continue; seen.add(k); uniq.push(s); if(uniq.length>=Math.min(P.maxSets,100))break; }
    // Estruturas secundárias: avalia hairpin/homo/heterodímero dos finalistas e
    // RE-RANQUEIA (conjuntos sem estrutura sobem). Roda só nos ≤100 finalistas (rápido).
    for (const s of uniq){
      const prims=[{name:'F3',seq:s.F3.seq},{name:'B3',seq:s.B3.seq},{name:'FIP',seq:s.FIP.seq},{name:'BIP',seq:s.BIP.seq}];
      if (s.LF) prims.push({name:'LF',seq:s.LF.seq});
      if (s.LB) prims.push({name:'LB',seq:s.LB.seq});
      const st=structureScan(prims, P.hpTh, P.dimerTh);
      s.penalty += st.penalty;
      s.struct={ hairpin:st.worstHairpin, dimer:st.worstDimer };
      if (st.warns.length) s.warns=(s.warns||[]).concat(st.warns);
    }
    uniq.sort((a,b)=>a.penalty-b.penalty);
    uniq.forEach((s,i)=>s.rank=i+1);
    return { sets:uniq, totalCandidates:sets.length, geom, params:P, targetLen:N, gc:Math.round(gcPct(S)), tReact:P.tReact };
  }

  // formata um aviso estruturado no idioma atual
  function fmtWarn(w) {
    const D=L();
    switch (w.t) {
      case 'gc':  return `${w.n} GC ${w.v}%`;
      case 'lc':  return `${w.n} LC ${w.v}%`;
      case 'len': return `${w.n} ${w.v} nt`;
      case 'dg':  return `${w.n} ${w.end}'ΔG ${w.v}`;
      case 'hairpin':   return `${w.n} ${D['warn.hairpin']} ΔG ${w.v}`;
      case 'selfdimer': return `${w.n} ${D['warn.selfdimer']} ΔG ${w.v}`;
      case 'dimer':     return `${w.n} ${D['warn.dimer']} ΔG ${w.v}`;
      default: return '';
    }
  }

  function renderResults(data) {
    lastData = data;
    resList.innerHTML = '';
    if (!data || !data.sets || data.sets.length === 0) {
      resEmpty.hidden = false; resList.hidden = true; return;
    }
    resEmpty.hidden = true; resList.hidden = false;

    const D = L();
    const P = data.params || {};
    const note = data.geom ? D['notes.geometry'](data.geom) : '';
    const header = document.createElement('div');
    header.className = 'lw-result-card';
    header.innerHTML = `
      <h4>${D['res.summary']}</h4>
      <div class="lw-kv"><div>${D['res.target']}</div><div>${data.targetLen} nt · GC ${data.gc}%</div></div>
      <div class="lw-kv"><div>${D['res.sets']}</div><div>${data.sets.length} ${D['res.setsOf'](data.totalCandidates)}</div></div>
      <div class="lw-kv"><div>${D['res.targetTm']}</div><div>F1c/B1c ${P.innerLo}–${P.innerHi} °C · F2/B2/F3/B3 ${P.outerLo}–${P.outerHi} °C</div></div>
      <div class="lw-kv"><div>${D['res.salts']}</div><div>Na⁺ ${P.na} mM · Mg²⁺ ${P.mg} mM · dNTP ${P.dntp} mM</div></div>
      <p class="lw-help">${note}</p>
      <p class="lw-help">${D['res.method']}</p>
    `;
    resList.appendChild(header);

    const copyBtn = (seq) => `<button class="lw-btn ghost lw-copy" data-seq="${seq}" title="${D['res.copy']}" style="padding:2px 8px; font-size:11px;">${D['res.copy']}</button>`;
    const row = (label, p) => p ? `<div class="lw-kv"><div>${label}</div><div>5'-${p.seq}-3' <span class="lw-suffix">(${p.len||p.seq.length} nt${p.tm!=null?`, Tm ${p.tm.toFixed(1)}°C`:''}${p.gc!=null?`, GC ${p.gc.toFixed(0)}%`:''})</span> ${copyBtn(p.seq)}</div></div>` : '';

    data.sets.forEach(set => {
      const card = document.createElement('div');
      card.className = 'lw-result-card';
      const warn = (set.warns&&set.warns.length)
        ? `<p class="lw-help" style="color:var(--warn);">⚠ ${set.warns.map(fmtWarn).join(' · ')}</p>`
        : `<p class="lw-help" style="color:var(--ok);">${D['res.allOk']}</p>`;
      card.innerHTML = `
        <h4>${D['res.set']} #${set.rank} — ${D['res.penalty']} ${set.penalty.toFixed(1)} · ${D['res.amplicon']} F3–B3 ${set.ampliconF3B3} nt</h4>
        <div class="lw-help" style="margin:0 0 4px;">${D['res.outer']}</div>
        ${row('F3', set.F3)}
        ${row('B3', set.B3)}
        <div class="lw-help" style="margin:8px 0 4px;">${D['res.inner']}</div>
        ${row('FIP', set.FIP)}
        ${row('↳ F1c', set.F1c)}
        ${row('↳ F2', set.F2)}
        ${row('BIP', set.BIP)}
        ${row('↳ B1c', set.B1c)}
        ${row('↳ B2', set.B2)}
        ${(set.LF||set.LB)?`<div class="lw-help" style="margin:8px 0 4px;">${D['res.loopHdr']}</div>`:''}
        ${row('LF', set.LF)}
        ${row('LB', set.LB)}
        <div class="lw-kv"><div>${D['res.positions']}</div><div>F3 ${set.coords.F3.join('–')} · F2 ${set.coords.F2.join('–')} · F1 ${set.coords.F1.join('–')} · B1 ${set.coords.B1.join('–')} · B2 ${set.coords.B2.join('–')} · B3 ${set.coords.B3.join('–')}</div></div>
        <div class="lw-kv"><div>${D['res.distances']}</div><div>F2–B2 ${set.ampliconF2B2} nt · F1–B1 ${set.dF1B1} nt</div></div>
        ${set.struct ? `<div class="lw-kv"><div>${D['res.struct'](data.tReact!=null?data.tReact:63)}</div><div>${D['res.worstHp']} ${set.struct.hairpin.toFixed(1)} · ${D['res.worstDi']} ${set.struct.dimer.toFixed(1)} kcal/mol</div></div>` : ''}
        ${warn}
      `;
      resList.appendChild(card);
    });

    resList.querySelectorAll('.lw-copy').forEach(b=>{
      b.addEventListener('click', ()=>{ if(navigator.clipboard) navigator.clipboard.writeText(b.dataset.seq); const t=b.textContent; b.textContent='✓'; setTimeout(()=>b.textContent=t,900); });
    });
  }

  btnExecutar.addEventListener('click', () => {
    const seq = sanitizeDNA(seqEl.value.trim());
    btnExecutar.textContent = L()['btn.generating'];
    btnExecutar.disabled = true;
    setTimeout(() => {
      try {
        const data = designLAMP(seq, readParams());
        if (data.errorKey) {
          alert(L()[data.errorKey]);
          renderResults({ sets: [] });
        } else if (!data.sets.length) {
          alert(L()['alert.noSets']);
          renderResults({ sets: [] });
        } else {
          renderResults(data);
          toTabResultados();
        }
      } catch (err) {
        console.error('LAMP design error:', err);
        alert(L()['alert.designErr'] + (err && err.message ? err.message : err));
      } finally {
        btnExecutar.textContent = L()['btn.run'];
        btnExecutar.disabled = false;
      }
    }, 30);
  });

  // ===================================================================
  // Especificidade — triagem in silico dos primers do Conjunto #1 contra
  // sequências de fundo: índice de k-mers + semente 3′ exata + extensão
  // sem gaps (modela anelamento iniciável; o 3′ deve parear para iniciar).
  // Privado: os primers não saem do navegador; só acessos públicos (se
  // usados) vão ao NCBI E-utilities para baixar as sequências de fundo.
  // ===================================================================
  const specAcc = qs('#spec-acc');
  const specBg = qs('#spec-bg');
  const specSeed = qs('#spec-seed');
  const specMaxmm = qs('#spec-maxmm');
  const btnFetch = qs('#btn-fetch-ncbi');
  const btnScreen = qs('#btn-screen');
  const specEmpty = qs('#spec-vazio');
  const specList = qs('#spec-list');
  const specOrg = qs('#spec-org');
  const specOrgList = qs('#spec-org-list');
  const btnListGenomes = qs('#btn-list-genomes');
  const specGenome = qs('#spec-genome');
  const btnAddGenome = qs('#btn-add-genome');
  let lastSpec = null; // {rank, primers, bgs, k, maxMM} — re-render ao trocar idioma

  function parseFastaMulti(text) {
    const out=[]; let cur=null;
    (text||'').split(/\r?\n/).forEach(line => {
      if (line.startsWith('>')) { if (cur) out.push(cur); cur={ name:(line.slice(1).trim() || ('seq'+(out.length+1))), seq:'' }; }
      else { const s=sanitizeDNA(line); if (cur) cur.seq+=s; else if (s) cur={ name:'seq1', seq:s }; }
    });
    if (cur) out.push(cur);
    return out.filter(r => r.seq.length>0).map(r => ({ name:r.name.slice(0,60), seq:r.seq }));
  }
  function buildKmerIndex(bg, k) {
    const idx=new Map();
    for (let i=0; i+k<=bg.length; i++){ const key=bg.substr(i,k); let a=idx.get(key); if(!a){a=[];idx.set(key,a);} a.push(i); }
    return idx;
  }
  function countMM(pat, bg, start, cap) {
    const Lp=pat.length; if (start<0 || start+Lp>bg.length) return cap+1;
    let mm=0; for (let i=0;i<Lp;i++){ if (bg[start+i]!==pat[i]){ if (++mm>cap) return mm; } }
    return mm;
  }
  // melhor alinhamento sem gaps de `primer` em `bg`, exigindo semente 3′ exata (k nt)
  function screenPrimerOnBg(primer, bg, idx, k, maxMM) {
    const Lp=primer.length; if (Lp<k) return null;
    let best=null;
    const consider=(ws, pat, strand)=>{ const mm=countMM(pat,bg,ws,maxMM); if (mm<=maxMM && (!best||mm<best.mm)) best={mm, start:ws, end:ws+Lp, strand}; };
    const occB=idx.get(primer.slice(Lp-k));      // primer ocorre em bg+ (3′ no fim do alinhamento)
    if (occB) for (const pos of occB) consider(pos-(Lp-k), primer, '+');
    const rc=revComp(primer);
    const occA=idx.get(rc.slice(0,k));           // revComp(primer) em bg+ (3′ no início)
    if (occA) for (const pos of occA) consider(pos, rc, '−');
    return best;
  }
  function primersFromSet(set) {
    const ps=[{name:'F3',seq:set.F3.seq},{name:'F2',seq:set.F2.seq},{name:'F1c',seq:set.F1c.seq},
              {name:'B1c',seq:set.B1c.seq},{name:'B2',seq:set.B2.seq},{name:'B3',seq:set.B3.seq}];
    if (set.LF) ps.push({name:'LF',seq:set.LF.seq});
    if (set.LB) ps.push({name:'LB',seq:set.LB.seq});
    return ps;
  }
  function renderSpec(spec) {
    lastSpec = spec; const D=L();
    specList.innerHTML='';
    if (!spec || !spec.bgs.length) { specEmpty.hidden=false; specList.hidden=true; return; }
    specEmpty.hidden=true; specList.hidden=false;
    const head=document.createElement('div'); head.className='lw-result-card';
    head.innerHTML=`<h4>${D['spec.setEval'](spec.rank)}</h4><p class="lw-help">${D['spec.note']}</p>`;
    specList.appendChild(head);
    const minLen=Math.min.apply(null, spec.primers.map(p=>p.seq.length));
    const k=Math.max(8, Math.min(spec.k, minLen));
    spec.bgs.forEach(bg => {
      const idx=buildKmerIndex(bg.seq, k);
      let hits=0; const rows=[];
      spec.primers.forEach(p => {
        const m=screenPrimerOnBg(p.seq, bg.seq, idx, k, spec.maxMM);
        if (m){ hits++; rows.push(`<p class="lw-help" style="margin:2px 0; color:var(--warn);">${D['spec.hit'](`${p.name} (5'-${p.seq}-3')`, m.mm, `${m.start+1}–${m.end}`, m.strand)}</p>`); }
        else rows.push(`<p class="lw-help" style="margin:2px 0; color:var(--ok);">${D['spec.clean'](p.name)}</p>`);
      });
      const card=document.createElement('div'); card.className='lw-result-card';
      card.innerHTML=`<h4>${D['spec.bgHdr'](bg.name, bg.seq.length)}</h4>`+
        `<p class="lw-help">${D['spec.summary'](hits, spec.primers.length)}</p>`+
        rows.join('')+
        (hits>=2 ? `<p class="lw-help" style="color:var(--warn);">${D['spec.multi']}</p>` : '');
      specList.appendChild(card);
    });
  }

  // E-utilities: baixa o FASTA (nuccore) de uma lista de acessos
  async function efetchFasta(ids) {
    const parts=[];
    for (const id of ids){
      const url=`https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${encodeURIComponent(id)}&rettype=fasta&retmode=text`;
      const r=await fetch(url); if (!r.ok) throw new Error(`${id}: HTTP ${r.status}`);
      const t=(await r.text()).trim(); if (!t.startsWith('>')) throw new Error(`${id}: ${t.slice(0,80)}`);
      parts.push(t);
    }
    return parts.join('\n');
  }
  function appendBg(text) { specBg.value=(specBg.value.trim() ? specBg.value.trim()+'\n' : '')+text; }

  if (btnFetch) btnFetch.addEventListener('click', async () => {
    const accs=(specAcc.value||'').split(/[\s,;]+/).filter(Boolean);
    if (!accs.length) { alert(L()['alert.specNoAcc']); return; }
    const orig=btnFetch.textContent; btnFetch.textContent=L()['btn.fetching']; btnFetch.disabled=true;
    try { appendBg(await efetchFasta(accs)); }
    catch (err) { alert(L()['alert.specFetchErr']+(err && err.message ? err.message : err)); }
    finally { btnFetch.textContent=orig; btnFetch.disabled=false; }
  });

  // Autocomplete de organismo — NCBI Datasets taxon_suggest -> <datalist> (terminologia NCBI)
  let orgTimer=null;
  function loadOrgSuggest(q) {
    q=(q||'').trim(); if (q.length<3 || !specOrgList) return;
    fetch('https://api.ncbi.nlm.nih.gov/datasets/v2alpha/taxonomy/taxon_suggest/'+encodeURIComponent(q)+'?tax_rank_filter=species')
      .then(r => r.ok ? r.json() : null)
      .then(j => {
        if (!j || !j.sci_name_and_ids) return;
        const seen=new Set(); specOrgList.innerHTML='';
        j.sci_name_and_ids.slice(0,12).forEach(s => {
          const n=s.sci_name; if (n && !seen.has(n)) { seen.add(n); const o=document.createElement('option'); o.value=n; specOrgList.appendChild(o); }
        });
      }).catch(()=>{});
  }
  if (specOrg) specOrg.addEventListener('input', () => { clearTimeout(orgTimer); orgTimer=setTimeout(() => loadOrgSuggest(specOrg.value), 300); });

  // Lista as sequências RefSeq genômicas do organismo (esearch + esummary, nuccore)
  if (btnListGenomes) btnListGenomes.addEventListener('click', async () => {
    const org=((specOrg && specOrg.value) || '').trim();
    if (!org) { alert(L()['alert.specNoOrg']); return; }
    const D=L(); const orig=btnListGenomes.textContent; btnListGenomes.textContent=D['btn.listing']; btnListGenomes.disabled=true;
    try {
      const term=encodeURIComponent('"'+org+'"[Organism] AND refseq[filter] AND biomol_genomic[PROP]');
      const r=await fetch('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nuccore&term='+term+'&retmax=20&retmode=json');
      const j=await r.json(); const ids=(j.esearchresult && j.esearchresult.idlist) || [];
      specGenome.innerHTML='';
      if (!ids.length) { const o=document.createElement('option'); o.value=''; o.textContent=D['spec.genomeNoneFound']; specGenome.appendChild(o); return; }
      const r2=await fetch('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=nuccore&id='+ids.join(',')+'&retmode=json');
      const j2=await r2.json();
      ids.forEach(id => {
        const it=j2.result && j2.result[id]; if (!it) return;
        const o=document.createElement('option'); o.value=it.accessionversion;
        o.textContent=D['spec.genomeOpt'](it.accessionversion, it.slen, (it.title||'').slice(0,60));
        specGenome.appendChild(o);
      });
    } catch (err) { alert(L()['alert.specGenomeErr']+(err && err.message ? err.message : err)); }
    finally { btnListGenomes.textContent=orig; btnListGenomes.disabled=false; }
  });

  // Adiciona o genoma/sequência escolhido ao campo de fundo (efetch)
  if (btnAddGenome) btnAddGenome.addEventListener('click', async () => {
    const acc=specGenome && specGenome.value;
    if (!acc) { alert(L()['alert.specNoGenome']); return; }
    const orig=btnAddGenome.textContent; btnAddGenome.textContent=L()['btn.fetching']; btnAddGenome.disabled=true;
    try { appendBg(await efetchFasta([acc])); }
    catch (err) { alert(L()['alert.specFetchErr']+(err && err.message ? err.message : err)); }
    finally { btnAddGenome.textContent=orig; btnAddGenome.disabled=false; }
  });

  if (btnScreen) btnScreen.addEventListener('click', () => {
    if (!lastData || !lastData.sets || !lastData.sets.length) { alert(L()['alert.specNoDesign']); return; }
    const bgs=parseFastaMulti(specBg.value);
    if (!bgs.length) { alert(L()['alert.specNoBg']); return; }
    const k=Math.max(8, Math.min(20, parseInt((specSeed && specSeed.value) || 13, 10) || 13));
    const maxMM=Math.max(0, Math.min(6, parseInt((specMaxmm && specMaxmm.value) || 2, 10) || 0));
    const set=lastData.sets[0];
    const primers=primersFromSet(set);
    const orig=btnScreen.textContent; btnScreen.textContent=L()['btn.screening']; btnScreen.disabled=true;
    setTimeout(() => {
      try { renderSpec({ rank:set.rank||1, primers, bgs, k, maxMM }); }
      catch (err) { console.error('spec error', err); alert(L()['alert.designErr']+(err && err.message ? err.message : err)); }
      finally { btnScreen.textContent=orig; btnScreen.disabled=false; }
    }, 30);
  });

  // Presets de exemplo (reatividade cruzada): só preenchem os campos;
  // o usuário ainda clica em "Buscar FASTA no NCBI" e depois "Triar".
  function fillSpecExample(acc, seed, maxmm) {
    if (specAcc) specAcc.value = acc;
    if (specSeed) specSeed.value = seed;
    if (specMaxmm) specMaxmm.value = maxmm;
  }
  (function () {
    const elA = document.getElementById('specEgAcentrale');
    if (elA) elA.addEventListener('click', () => fillSpecExample('CP001759.1', 13, 2));
    const elB = document.getElementById('specEgBbigemina');
    if (elB) elB.addEventListener('click', () => fillSpecExample('KP710228.1', 13, 2));
  })();

  // ===================================================================
  // Varredura de Mg²⁺ → heatmap (ADITIVO). Usa o conjunto desenhado atual
  // (mesmo que renderSpec/btnScreen: lastData.sets[0]) e o fundo carregado
  // (specBg). Para cada primer e cada MgSO₄ total, calcula a margem de Tm
  // on−off (tmDuplexMg, caminho separado) e sugere o [Mg²⁺] por maximin.
  // ===================================================================
  const specDntp = qs('#spec-dntp');
  const specMon = qs('#spec-mon');
  const specTrxn = qs('#spec-trxn');
  const specMgmin = qs('#spec-mgmin');
  const specMgmax = qs('#spec-mgmax');
  const specMgstep = qs('#spec-mgstep');
  const btnMgSweep = qs('#btn-mgsweep');
  const specMgHeatmap = qs('#specMgHeatmap');
  const specMgSuggest = qs('#specMgSuggest');
  let lastMg = null; // último resultado da varredura (para re-render ao trocar idioma)

  // faixa de Tm on-target válida por papel do primer (janela LAMP)
  function lampTmWindow(name){
    if (/^(F3|B3|F2|B2)$/.test(name)) return [59,61];
    if (/^(F1c|B1c|LF|LB)$/.test(name)) return [64,66];
    return [59,66]; // relaxado
  }
  function fmtMg(v){ return (Math.round(v*100)/100); }

  // Núcleo: varre Mg; por célula calcula S(i,j)=Tm_on−Tm_off (discriminação
  // intrínseca, ~invariante ao Mg) E safety(i,j)=T_rxn−Tm_off (margem de
  // segurança do off-target, DEPENDENTE do Mg). "functional" = Tm_on na janela
  // LAMP. Escolhe j* pela menor coluna de Mg em que todos são functional E safe.
  function computeMgSweep(primers, bgs, opt){
    const mgMin=opt.mgMin, mgMax=opt.mgMax, mgStep=opt.mgStep>0?opt.mgStep:0.5;
    const cols=[]; for (let mg=mgMin; mg<=mgMax+1e-9; mg+=mgStep) cols.push(Math.round(mg*1000)/1000);
    // off-target alinhado uma vez por primer (independe do Mg)
    const offs = primers.map(p => bestOffWindow(p.seq, bgs));
    const rows = primers.map((p,i)=>{
      const win=lampTmWindow(p.name);
      const cells = cols.map(mg=>{
        const o={ mg, dntp:opt.dntp, mon:opt.mon, oligoNM:TM_OLIGO_NM };
        const on = tmDuplexMg(p.seq, wcBottom(p.seq), o);
        let off={tm:NaN, flag:true};
        if (offs[i]) off = tmDuplexMg(p.seq, offs[i].alignedTargetBottom, o);
        const S = (isNaN(on.tm)||isNaN(off.tm)) ? NaN : (on.tm - off.tm);
        // margem de segurança do off-target (Mg-dependente): T_rxn − Tm_off.
        // Sem off-target (Tm_off=NaN) ⇒ infinitamente seguro.
        const safety = isNaN(off.tm) ? Infinity : (opt.tRxn - off.tm);
        const functional = !isNaN(on.tm) && inRange(on.tm, win[0], win[1]);
        const safe = safety > 0;
        // "valid" mantido p/ retrocompat: functional E safe (igual à regra antiga)
        const valid = functional && safe;
        const freeMg = freeMgM(mg/1000, opt.dntp/1000);
        return { mg, freeMg, tmOn:on.tm, tmOff:off.tm, S, safety, functional, safe, valid };
      });
      return { name:p.name, seq:p.seq, off:offs[i], cells };
    });
    // por coluna: minSafety=min_i safety; nPass=#(functional E safe); allPass
    const colStats = cols.map((mg,j)=>{
      let minSafety=Infinity, minS=Infinity, allPass=true, nPass=0;
      rows.forEach(r=>{
        const c=r.cells[j];
        if (!isNaN(c.S)) minS=Math.min(minS,c.S);
        if (isFinite(c.safety)) minSafety=Math.min(minSafety,c.safety);
        if (c.functional && c.safe) nPass++; else allPass=false;
      });
      if (!isFinite(minS)) minS=NaN;
      if (!isFinite(minSafety)) minSafety=NaN;
      return { mg, freeMg:freeMgM(mg/1000, opt.dntp/1000), minS, minSafety, feasible:allPass, nPass, nValid:nPass };
    });
    // j* = MENOR coluna de Mg com todos functional E safe; senão a coluna que
    // maximiza nPass, menor Mg em empate.
    let jStar=-1;
    for (let j=0;j<colStats.length;j++){ if (colStats[j].feasible){ jStar=j; break; } }
    if (jStar<0){
      let bestN=-1; colStats.forEach((c,j)=>{ if (c.nPass>bestN){ bestN=c.nPass; jStar=j; } });
    }
    return { cols, rows, colStats, jStar, feasible: jStar>=0 && colStats[jStar].feasible, opt };
  }

  // diverging centrado em 0: vermelho(safety≤0)→branco(0)→verde(safety>0).
  // Cinza p/ célula não-funcional (Tm_on fora da janela LAMP) ou Tm indefinida.
  function colorForSafety(safety, functional){
    if (isNaN(safety)) return '#9aa0a6';
    if (!functional) return '#c3c7cc';
    const span=10; // °C de escala total (clamp em ±10)
    const t=Math.max(-1, Math.min(1, safety/span));
    let r,g,b;
    if (t<=0){ const k=t+1; r=Math.round(198+(255-198)*k); g=Math.round(40+(255-40)*k); b=Math.round(40+(255-40)*k); }
    else { const k=1-t; r=Math.round(40+(255-40)*k); g=Math.round(167+(255-167)*k); b=Math.round(70+(255-70)*k); }
    return `rgb(${r},${g},${b})`;
  }

  function drawMgHeatmap(res){
    if (!specMgHeatmap || !specMgHeatmap.getContext) return;
    const ctx=specMgHeatmap.getContext('2d');
    const nR=res.rows.length, nC=res.cols.length;
    const padL=58, padT=34, cellW=Math.max(26, Math.min(54, 720/Math.max(1,nC))), cellH=26;
    const W=padL+nC*cellW+12, H=padT+nR*cellH+18;
    specMgHeatmap.width=W; specMgHeatmap.height=H;
    ctx.clearRect(0,0,W,H);
    ctx.font='11px system-ui,Arial,sans-serif'; ctx.textBaseline='middle';
    // cabeçalho das colunas: totMg / ~freeMg
    res.cols.forEach((mg,j)=>{
      const x=padL+j*cellW; const cs=res.colStats[j];
      ctx.save(); ctx.translate(x+cellW/2, padT-6); ctx.rotate(-Math.PI/4);
      ctx.fillStyle='#333'; ctx.textAlign='right';
      ctx.fillText(`${fmtMg(mg)}/${fmtMg(cs.freeMg*1000)}`, 0, 0);
      ctx.restore();
    });
    // linhas
    res.rows.forEach((r,i)=>{
      const y=padT+i*cellH;
      ctx.fillStyle='#222'; ctx.textAlign='right'; ctx.fillText(r.name, padL-8, y+cellH/2);
      r.cells.forEach((c,j)=>{
        const x=padL+j*cellW;
        // cor = margem de segurança do off-target (T_rxn − Tm_off), diverging em 0
        ctx.fillStyle=colorForSafety(c.safety, c.functional);
        ctx.fillRect(x+1, y+1, cellW-2, cellH-2);
        // anotação = Tm_off absoluto (°C); cinza se não-funcional
        if (!isNaN(c.tmOff)){
          ctx.fillStyle = c.functional ? '#102010' : '#555';
          ctx.textAlign='center';
          ctx.fillText(c.tmOff.toFixed(1), x+cellW/2, y+cellH/2);
        }
      });
    });
    // caixa de destaque na coluna j*
    if (res.jStar>=0){
      const x=padL+res.jStar*cellW;
      ctx.strokeStyle='#1558d6'; ctx.lineWidth=2.5;
      ctx.strokeRect(x+0.5, padT+0.5, cellW-1, nR*cellH-1);
    }
  }

  function renderMgSuggest(res){
    if (!specMgSuggest) return; const D=L();
    if (!res){ specMgSuggest.innerHTML=''; return; }
    const j=res.jStar; const cs = (j>=0) ? res.colStats[j] : null;
    let html='';
    if (cs){
      const tot=fmtMg(cs.mg), free=fmtMg(cs.freeMg*1000), minSafety=isNaN(cs.minSafety)?'—':cs.minSafety.toFixed(1);
      html += `<p class="lw-help" style="font-weight:600;">${res.feasible ? D['spec.mg.suggest'](tot, free, minSafety) : D['spec.mg.noFeas'](tot, free)}</p>`;
    } else {
      html += `<p class="lw-help">${D['spec.mg.none']}</p>`;
    }
    // tabela por primer em j*: #mm | S=Tm_on−Tm_off | Tm_on | Tm_off | safety | status
    if (j>=0){
      const head=`<tr><th style="text-align:left;">${D['spec.mg.primer']}</th><th>#mm</th><th>S (Δ)</th><th>Tm on</th><th>Tm off</th><th>${D['spec.mg.safety']}</th><th>${D['spec.mg.status']}</th></tr>`;
      const body=res.rows.map(r=>{
        const c=r.cells[j];
        const onS=isNaN(c.tmOn)?'—':c.tmOn.toFixed(1), offS=isNaN(c.tmOff)?'—':c.tmOff.toFixed(1);
        const sS=isNaN(c.S)?'—':((c.S>=0?'+':'')+c.S.toFixed(1));
        const safeS=isFinite(c.safety)?((c.safety>=0?'+':'')+c.safety.toFixed(1)):'∞';
        const mmS = r.off ? String(r.off.nMismatch) : '—';
        // status: ✓ safe & functional; ✗ caso contrário, com a falha apontada
        let mark;
        if (c.functional && c.safe) mark='<span style="color:var(--ok);">✓</span>';
        else { const why=[]; if (!c.functional) why.push(D['spec.mg.notFunc']); if (!c.safe) why.push(D['spec.mg.notSafe']);
          mark=`<span style="color:var(--warn);">✗ ${why.join(', ')}</span>`; }
        return `<tr><td style="text-align:left;">${r.name}</td><td>${mmS}</td><td>${sS}</td><td>${onS}</td><td>${offS}</td><td>${safeS}</td><td>${mark}</td></tr>`;
      }).join('');
      html += `<table class="lw-mg-table" style="width:100%; border-collapse:collapse; font-size:12px;"><thead>${head}</thead><tbody>${body}</tbody></table>`;
    }
    if (!res.feasible){
      const off=res.rows.filter(r=>{ const c=r.cells[Math.max(0,j)]; return !(c.functional && c.safe); }).map(r=>r.name);
      if (off.length) html += `<p class="lw-help" style="color:var(--warn);">${D['spec.mg.offenders'](off.join(', '))}</p>`;
    }
    html += `<p class="lw-help" style="font-style:italic;">${D['spec.mg.note']}</p>`;
    html += `<p class="lw-help" style="font-style:italic;">${D['spec.mg.note2']}</p>`;
    specMgSuggest.innerHTML=html;
  }

  function buildMgHeatmap(){
    const D=L();
    if (!lastData || !lastData.sets || !lastData.sets.length){ alert(D['alert.specNoDesign']); return; }
    const bgs=parseFastaMulti(specBg ? specBg.value : '');
    if (!bgs.length){ alert(D['alert.specNoBg']); return; }
    const num=(el,d)=>{ const v=parseFloat(el&&el.value); return isNaN(v)?d:v; };
    const opt={
      dntp:num(specDntp,1.6), mon:num(specMon,50), tRxn:num(specTrxn,63),
      mgMin:num(specMgmin,2), mgMax:num(specMgmax,10), mgStep:num(specMgstep,0.5)
    };
    const set=lastData.sets[0];
    const primers=primersFromSet(set);
    const res=computeMgSweep(primers, bgs, opt);
    lastMg=res;
    drawMgHeatmap(res);
    renderMgSuggest(res);
  }

  if (btnMgSweep) btnMgSweep.addEventListener('click', () => {
    try { buildMgHeatmap(); }
    catch (err){ console.error('mg sweep error', err); alert(L()['alert.designErr']+(err && err.message ? err.message : err)); }
  });

  // presets de dNTP total (baixo 1.6 mM / alto 5.6 mM)
  (function(){
    const a=document.getElementById('specDntpPE'); if (a && specDntp) a.addEventListener('click', ()=>{ specDntp.value='1.6'; });
    const b=document.getElementById('specDntpNEB'); if (b && specDntp) b.addEventListener('click', ()=>{ specDntp.value='5.6'; });
  })();

  // ===================================================================
  // i18n — aplicação no DOM + botão de idioma
  // ===================================================================
  function applyLang(lang) {
    LANG = (lang === 'en') ? 'en' : 'pt';
    const D = L();
    qsa('[data-i18n]').forEach(el => { const k=el.getAttribute('data-i18n'); if (typeof D[k]==='string') el.textContent = D[k]; });
    qsa('[data-i18n-html]').forEach(el => { const k=el.getAttribute('data-i18n-html'); if (typeof D[k]==='string') el.innerHTML = D[k]; });
    qsa('[data-i18n-ph]').forEach(el => { const k=el.getAttribute('data-i18n-ph'); if (typeof D[k]==='string') el.setAttribute('placeholder', D[k]); });
    const tg = qs('#langToggle');
    if (tg) { tg.textContent = (LANG==='en') ? 'Versão em Português' : 'English Version'; tg.title = tg.textContent; }
    document.documentElement.lang = (LANG==='en') ? 'en' : 'pt-BR';
    // botões com estado padrão acompanham o idioma
    if (btnExecutar && !btnExecutar.disabled) btnExecutar.textContent = D['btn.run'];
    if (btnAplicar) btnAplicar.textContent = D['btn.apply'];
    // re-renderiza resultados em cache no novo idioma
    if (lastData) renderResults(lastData);
    if (lastSpec) renderSpec(lastSpec);
    if (lastMg) renderMgSuggest(lastMg);
  }

  // O botão de idioma é um LINK NATIVO entre LAMPrime.html (PT) e LAMPrime_en.html (EN);
  // não há listener de toggle. O idioma de cada página vem do <html lang>.
  applyLang(LANG);
})();
