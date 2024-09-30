import glob
import pandas as pd

class GenUInE:
    def __init__(self, path, outputName="ANALYSIS", genome="hg19.bed", window=10000):
        self.path = path
        self.outputName = outputName
        self.genome = genome
        self.window = window

    @staticmethod
    def range_sequence(start, stop, step):
        """Create a list of multiple ranges based on a window"""
        result = list(zip(range(start, stop, step), range(start + step - 1, stop, step)))
        if (stop - start) % step != 0:
            last_fst_elem = result[-1][-1] if result else start
            result.append((last_fst_elem + 1, stop))
        else:
            result = result[:-1]
            last_fst_elem = result[-1][-1] if result else start
            result.append((last_fst_elem + 1, stop))
        return result

    @staticmethod
    def range_subset(range1, range2):
        """Check if range1 is a subset of range2."""
        if range1 and range2:
            return range1.start in range2 and (range1.stop - 1) in range2
        return False

    @staticmethod
    def range_extreme1(range1, range2):
        return range1.start in range2

    @staticmethod
    def range_extreme2(range1, range2):
        return range1.stop in range2

    @staticmethod
    def apply_stat(df, num):
        dfcopy = df.copy()
        for val in range(1, num + 1):
            pwin = dfcopy.iloc[:, val].value_counts().get(1, 0) / len(df)
            dfcopy.iloc[:, val] = dfcopy.iloc[:, val].replace(1, pwin)
            dfcopy.iloc[:, val] = dfcopy.iloc[:, val].replace(0, 1)
        df["comb_prob_value"] = dfcopy.iloc[:, 1:num + 1].prod(axis=1).astype(float)
        alpha = 0.5  # weight for combined probability
        beta = 0.5   # weight for summation
        df["score"] = alpha * (-10*np.log(df["comb_prob_value"])) + beta * (df["Summation"] / num)
        return df

    def matrix_gen(self):
        """Generate binary matrix every defined window"""

        print("Create reference intervals...\n")

        ref_raw = pd.read_csv(self.genome, sep="\t", names=["chr", "start", "end"])
        ref_coord = {}
        for chrom, start, end in zip(ref_raw["chr"], ref_raw["start"], ref_raw["end"]):
            ref_coord[chrom] = GenUInE.range_sequence(start, end, self.window)

        # Convert the ranges into sets for faster lookup
        for e, v in ref_coord.items():
            ref_coord[e] = [range(elem[0], elem[1]) for elem in v]

        df_ref = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in ref_coord.items()]))
        for col in df_ref:
            df_ref[col] = f'{col}:' + df_ref[col].astype(str)

        df_ref = pd.concat([df_ref, df_ref.T.stack().reset_index(name='Genome')['Genome']], axis=1)
        df_ref = df_ref[["Genome"]]
        df_ref = df_ref[~df_ref["Genome"].str.contains("nan")]

        print("Reference uploaded!\nCreate matrix with selected .bed\n")

        n_bed = 0
        for file in glob.glob(f"{self.path}*.bed"):
            fileName = file.split("/")[-1].split(".")[0]
            print(fileName)
            n_bed += 1
            df_bed = pd.read_csv(file, sep="\t", names=["chr", "start", "end"])
            df_dict = df_bed.groupby('chr').apply(lambda x: list(zip(x['start'], x['end']))).to_dict()

            # Pre-calculate ranges for df_dict to avoid recalculation in loops
            range_cache = {}
            for e, v in df_dict.items():
                df_dict[e] = [range(elem[0], elem[1]) for elem in v]
                for i in df_dict[e]:
                    if len(i) > self.window:
                        range_cache[i] = GenUInE.range_sequence(i.start, i.stop, self.window)

            data = set()

            # Iterare su tutti i cromosomi comuni tra ref_coord e df_dict
            common_chromosomes = set(ref_coord.keys()).intersection(df_dict.keys())
            for chr in common_chromosomes:
                ref_ranges = ref_coord[chr]
                bed_ranges = df_dict[chr]

                for ref_range in ref_ranges:
                    found_match = False

                    for bed_range in bed_ranges:
                        if len(bed_range) > self.window:
                            rangelist = range_cache.get(bed_range, [])
                           # print(rangelist)
                            if any(GenUInE.range_subset(range(r[0], r[1]), ref_range) or
                                   GenUInE.range_extreme1(range(r[0], r[1]), ref_range) or
                                   GenUInE.range_extreme2(range(r[0], r[1]), ref_range) for r in rangelist):
                                data.add(f"{chr}:{ref_range}")
                                print(f"{chr}:{ref_range}")
                                found_match = True
                                break
                        else:
                            if GenUInE.range_subset(bed_range, ref_range) or \
                               GenUInE.range_extreme1(bed_range, ref_range) or \
                               GenUInE.range_extreme2(bed_range, ref_range):
                                data.add(f"{chr}:{ref_range}")
                                print(f"{chr}:{ref_range}")
                                found_match = True
                                break

                    if found_match:
                        continue

            df_mat = pd.DataFrame(list(data), columns=["Genome"])
            df_mat[fileName] = 1
            df_ref = df_ref.merge(df_mat, on="Genome", how="left").fillna(0)

        df_ref['Summation'] = df_ref.iloc[:, 1:].sum(axis=1)
        df_ref = GenUInE.apply_stat(df_ref, n_bed)
        df_ref.to_csv(f"ENRICHED_MATRIX_{self.outputName}.tsv", sep='\t', index=False)
        print("Enriched matrix created!")