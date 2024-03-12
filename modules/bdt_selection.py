class DQCDBDTSelectionRDFProducer():
    def __init__(self, *args, **kwargs):
        self.bdt_name = kwargs.pop("bdt_name", "bdt_scenarioA")
        self.bdt_cut_value = kwargs.pop("bdt_cut_value", 0.)

    def run(self, df):
        df = df.Filter(f"{self.bdt_name} > {self.bdt_cut_value}", f"bdt > {self.bdt_cut_value}")
        return df, []


def DQCDBDTSelectionRDF(*args, **kwargs):
    return lambda: DQCDBDTSelectionRDFProducer(*args, **kwargs)