``` mermaid
flowchart TD
    subs["Subjects \n 11,878"]
    %% geno[Unrelated Genotypes \n 10,635 (-10%)]
    %% unrel[Unrelated Genotypes \n 7,234 (-32%)]
    unrel["Unrelated Genotypes \n 7,234 (-39%)"]

    fMRI["Gordon \n 5,709 (-52%)"]

    herit["SNP Heritabilty \n 4,850 (-15%)"]

    subs --> unrel --> herit
    subs --> fMRI --> herit

