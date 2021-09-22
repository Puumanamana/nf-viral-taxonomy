import argparse
import pandas as pd


RANKS = ['superkingdom', 'realm', 'kingdom', 'phylum', 'class',
         'order', 'family', 'subfamily', 'genus', 'species', 'strain']

def parse_args():

    parser = argparse.ArgumentParser()
    parser.add_argument('--hits', type=str)
    parser.add_argument('--mapping', type=str)
    parser.add_argument('--outdir', type=str, default='.')
    parser.add_argument('--vote', type=float, default=0.5)
    parser.add_argument('--min-hits', type=float, default=2)
    args = parser.parse_args()

    return args

def clean_rank(entry):
    """
    match <names> entries with <RANKS> and fill the blanks in <names> (="clade" or "no rank") 
    if the number of blanks matches the number of missing ranks in <RANKS>
    # Account for blank at the end?
    """

    if isinstance(entry, float):
        return entry

    fixed_ranks = []
    diff = 0
    last_known_index = -1
    for i, x in enumerate(entry):
        if x not in RANKS:
            diff += 1
            continue

        if diff > 0 and i-last_known_index-1 == diff:
            fixed_ranks += RANKS[last_known_index+1:i]

        fixed_ranks.append(x)
        diff = 0
        last_known_index = i

    return fixed_ranks

def main():
    args = parse_args()

    contig_hits = pd.read_csv(args.hits, sep='\t', usecols=['scaffold', 'viral_id'])
    contig_hits.viral_id = contig_hits.viral_id.str.split('.').str[0]
    
    hits_mapping = pd.read_csv(args.mapping, sep='\t', header=None, names=['viral_id', 'taxid', 'lineage', 'ranks'])

    hits_mapping = hits_mapping.dropna()
    hits_mapping.ranks = (hits_mapping.ranks
                  .str.lower()
                  .str.replace(r'genus;no rank$', 'genus;species', regex=True)
                  .str.replace(r'species;no rank$', 'species;strain', regex=True)
                  .str.split(';')
                  .map(clean_rank))
    hits_mapping.lineage = hits_mapping.lineage.str.split(';')
    hits_mapping.lineage = hits_mapping[['ranks', 'lineage']].apply(lambda x: ['|'.join(y) for y in zip(*x)], axis=1)

    lineages = (
        contig_hits.dropna()
        .merge(hits_mapping, on='viral_id')
        .set_index('scaffold')
        .lineage
        .explode()
        .str.split('|', expand=True)
    )
    
    lineages.columns = ['rank', 'value']

    # Remove contigs with only 1 protein called
    group_sizes = lineages[lineages['rank']=='superkingdom'].groupby(level='scaffold').size()
    # lineages = lineages[group_sizes.loc[lineages.index]>=args.min_hits].reset_index().set_index('rank')
    lineages = lineages.reset_index().set_index('rank')

    assignment = dict()
    for rank in RANKS:
        # Take the first count in value_counts (largest count is always first)
        votes = lineages.loc[rank].value_counts().rename('freq').reset_index()
        votes['freq'] /= group_sizes.loc[votes.scaffold].values
        # majority_vote = votes.groupby('qseqid').first()
        majority_vote = votes.groupby('scaffold').apply(lambda x: x.loc[x.freq.idxmax()])
        majority_vote = majority_vote[majority_vote.freq > args.vote]
        assignment[rank] = majority_vote['value']

    assignment = pd.DataFrame(assignment).reindex(contig_hits.scaffold.unique())
    print(assignment.count())

    assignment.to_csv(f'{args.outdir}/taxonomy.csv')

if __name__ == '__main__':
    main()
