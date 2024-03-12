import concurrent.futures


class GeneAnnotationsSearch:
    # ... (existing code remains the same)

    def get_gene_annotations_bulk(self, regions: List[str], chunk_size: int = 1000,
                                  num_threads: int = os.cpu_count()) -> List[List[str]]:
        """
        :param regions: list of regions to search for, e.g. ["1:1000-1000000", "2:2000-2000000"]
        :param chunk_size: number of regions to process in each chunk (default: 1000)
        :param num_threads: number of threads to use for parallel processing (default: number of CPU cores)
        :return: list of lists of annotations, corresponding to each region
        """
        annotations = [[] for _ in regions]

        def process_chunk(chunk_start: int, chunk_end: int) -> List[List[str]]:
            chunk_regions = regions[chunk_start:chunk_end]
            chunk_annotations = self._get_gene_annotations_chunk(chunk_regions)
            return chunk_annotations

        chunk_starts = list(range(0, len(regions), chunk_size))
        chunk_ends = chunk_starts[1:] + [len(regions)]

        with concurrent.futures.ThreadPoolExecutor(max_workers=num_threads) as executor:
            futures = [
                executor.submit(process_chunk, start, end)
                for start, end in zip(chunk_starts, chunk_ends)
            ]
            results = [future.result() for future in concurrent.futures.as_completed(futures)]

        for chunk_annotations in results:
            for i, annotation in enumerate(chunk_annotations):
                annotations[chunk_starts[i]] = annotation
                chunk_starts[i] += 1

        return annotations

    def _get_gene_annotations_chunk(self, regions: List[str]) -> List[List[str]]:
        """
        :param regions: list of regions to search for, e.g. ["1:1000-1000000", "2:2000-2000000"]
        :return: list of lists of annotations, corresponding to each region
        """
        annotations = [[] for _ in regions]

        query = " OR ".join(
            f"(seqid = '{region.split(':')[0]}' AND start <= {region.split(':')[1].split('-')[1]} AND end >= {region.split(':')[1].split('-')[0]})"
            for region in regions)

        results = self.db.execute(f"""
            SELECT seqid, start, end, featuretype
            FROM features
            WHERE {query}
        """)

        for result in results:
            for i, region in enumerate(regions):
                if result.seqid == region.split(':')[0] and result.start <= int(
                        region.split(':')[1].split('-')[1]) and result.end >= int(region.split(':')[1].split('-')[0]):
                    if result.featuretype not in annotations[i]:
                        annotations[i].append(result.featuretype)

        for i in range(len(annotations)):
            annotations[i].sort()

        return annotations