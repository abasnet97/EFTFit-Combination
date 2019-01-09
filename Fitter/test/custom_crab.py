def custom_crab(config):
    print '>> Customising the crab config'
    from CRABClient.UserUtilities import getUsernameFromSiteDB
    config.Site.storageSite = 'T3_US_NotreDame'
    config.Site.blacklist = ['T3_IT_Bologna', 'T3_US_UMiss', 'T2_ES_IFCA', 'T2_TR_METU', 'T2_CH_CSCS', 'T3_US_Baylor']
    config.Data.outLFNDirBase = '/store/user/%s/EFT/' % (getUsernameFromSiteDB())
    #config.Data.outLFNDirBase = '/store/user/alefeld/EFT/'